import subprocess
from pathlib import Path
import pandas as pd
import numpy as np
from Bio import SeqIO, Seq, SeqRecord
import os
from glob import glob

def run_sys_command(command):
    """
    Run a shell command using subprocess.

    Usage:
    run_sys_command(command)

    Inputs:
    command: str, shell command
    
    Returns:
    None
    """
    print(f'Running command: {command}')
    result = subprocess.run(
        command, 
        shell=True, 
        stdout=subprocess.DEVNULL,  # Suppress stdout
        stderr=subprocess.PIPE      # Capture stderr, if needed
    )
    
    # Check the return code for success or failure
    if result.returncode == 0:
        print("Command executed successfully.")
    else:
        print("Command failed. Error:", result.stderr.decode())
    return
    
def get_fastq_read_count(fastq_path):
    """
    Count the number of reads a fastq.

    Usage:
    get_fastq_read_count(fastq_path)

    Inputs:
    fastq_path: str or path, path to fastq file
    
    Returns:
    read_count: int, number of reads in fastq file
    """
    # Run the wc command and capture output
    command = f"wc -l {fastq_path} | awk '{{print $1/4}}'"
    result = subprocess.run(command, shell=True, capture_output=True, text=True)
    
    # Extract the line count and compute the number of reads
    read_count = int(result.stdout)  # Extract first number (line count)
    
    return read_count

def parse_base_string(basestring, qualstring):
    """
    Convert the pileup basestring and qualstring information into friendlier lists.

    Usage:
    parse_base_string(basestring, qualstring)

    Inputs:
    basestring: str, pileup read string, ie. '...A..T.'
    qualstring: str, qscore string from pileup
    
    Returns:
    dataframe with position index and columns with the base identity, qscore, and sense
    """
    phred_dict = {chr(score + 33): score for score in range(0, 94)}
    
    baselist = []
    quallist = []
    senselist = []
    bind = 0
    qind = 0

    while bind < len(basestring):
        char = basestring[bind]
        # For simple aligned bases (matches and mismatches)
        if char in '.ATCG':
            sense = '+'
            baselist.append(char)
            quallist.append(phred_dict[qualstring[qind]])
            senselist.append(sense)
            bind += 1
            qind += 1
            continue
        elif char in ',atcg':
            sense = '-'
            baselist.append(char)
            quallist.append(phred_dict[qualstring[qind]])
            senselist.append(sense)
            bind += 1
            qind += 1
            continue            
        elif char == '$':  
            # End-of-read marker; usually no quality associated, so you might skip it
            # baselist.append('$')
            bind += 1
            continue
        elif char == '*':
            baselist.append('*')
            quallist.append(pd.NA)
            # You may choose to assign a strand that reflects a deletion – 
            # for example, if you want to mark it as 'del' or copy the last known strand.
            senselist.append(pd.NA)
            bind += 1
            continue

        elif char in "+-":
            # Handle insertions/deletions
            indel_symbol = char
            indel = [indel_symbol]
            bind += 1
            # Get all digits that represent the indel length
            num_str = ""
            while bind < len(basestring) and basestring[bind].isdigit():
                num_str += basestring[bind]
                bind += 1
            if num_str == "":
                raise ValueError("Expected indel length but found nothing.")
            indel.append(num_str)
            try:
                indel_length = int(num_str)
            except ValueError:
                raise ValueError(f"Invalid indel length: {num_str}")
            # Append the inserted/deleted sequence
            for _ in range(indel_length):
                if bind < len(basestring):
                    indel.append(basestring[bind])
                    bind += 1
                else:
                    raise ValueError("Indel length exceeds available bases.")
            if indel[-1] in '.ATCG':
                sense = '+'
            elif indel[-1] in ',atcg':
                sense = '-'
            else:
                raise ValueError(f"Problem detecting sense at {bind}")
            senselist.append(sense)   
            baselist.append(''.join(indel))
            quallist.append(pd.NA)
            continue
        elif char == '^':
            # Start-of-read marker plus mapping quality; we can either skip or record it.
            # Here, we'll record it as a token without a quality value.
            token = basestring[bind:bind+2]
            # baselist.append(token)
            # quallist.append(pd.NA)
            bind += 2
            continue
        else:
            # Unexpected character; you might decide to raise an error or skip.
            print(f"Warning: Unexpected character '{char}' in pileup string.")
            bind += 1
    
    return pd.DataFrame.from_dict({'base':baselist, 'qscore': quallist, 'sense': senselist})

def paired_to_consensus(r1_fastq, r2_fastq, consensus_fastq):
    """
    Generate a consensus fastq from paired end read Illumina fastq files.

    Usage:
    paired_to_consensus(r1_fastq, r2_fastq, consensus_fastq)

    Description:
    Runs the usearch tool for merging paired end reads into a single consensus fastq

    Inputs:
    r1_fastq: str, path to R1 fastq
    r2_fastq: str, path to R2 fastq
    consensus_fastq: str, path to output consensus file
    
    Returns:
    None
    """
    usearch_call = f"usearch -fastq_mergepairs {r1_fastq} -reverse {r2_fastq} -fastqout {consensus_fastq} -relabel @"
    run_sys_command(usearch_call)
    return

def extract_rxn_fastq(fastq_file, for_primer0, rep_primer0, out_fastq, spacer_for=True):
    """
    Extracts a fastq for a particular primer pair to output the subset of reads matching those primers. 

    Usage:
    extract_rxn_fastq(fastq_file, for_primer0, rep_primer0, out_fastq, spacer_for=True)

    Inputs:
    fastq_file: str or path, path to consensus fastq file with mixed reads
    for_primer0: str or biopython seq, sequence of forward primer for amplicon
    rep_primer0: str or biopython seq, sequence of reverse primer for amplicon
    out_fastq: str or path, name of output fastq file for filtered reads matching the primers.
    spacer_for: boolean, True if spacer and forward primer are on the same strand, False if on different strands
    
    Returns:
    None
    """

    # get the forward primer into biopython sequence form, generate its reverse complement for reverse reads.
    # also get the reverse complement of the reverse primer for the PCR with the replicate barcode.
    for_primer = Seq.Seq(for_primer0)
    rev_primer = for_primer.reverse_complement()
    rep_primer = Seq.Seq(rep_primer0).reverse_complement()

    # generate paths to temporary fastq files
    for_temp_fastq = Path(fastq_file).with_name('for_temp.fastq')
    rev_temp_fastq = Path(fastq_file).with_name('rev_temp.fastq')

    # use cutadapt in order to extract the reads matching for_primer in the forward sense and in the reverse sense
    for_cutadapt_call = f"cutadapt -g ^{for_primer} --discard-untrimmed -o {for_temp_fastq} {fastq_file}"
    rev_cutadapt_call = f"cutadapt -a {rev_primer}$ --discard-untrimmed -o {rev_temp_fastq} {fastq_file}"

    run_sys_command(for_cutadapt_call)
    run_sys_command(rev_cutadapt_call)

    print(f"Found {get_fastq_read_count(for_temp_fastq)} forward reads and {get_fastq_read_count(rev_temp_fastq)} reverse reads.")

    # take the fastq matching the reverse sense and flip it to the forward sense.
    rev_temp_revcomp_fastq = rev_temp_fastq.with_name('rev_temp_revcomp.fastq')
    run_sys_command(f"seqtk seq -r {rev_temp_fastq} > {rev_temp_revcomp_fastq}")    

    # append the flipped sense reads and remove temporary files
    os.system(f"cat {str(rev_temp_revcomp_fastq)} >> {for_temp_fastq}")
    os.remove(rev_temp_fastq)
    os.remove(rev_temp_revcomp_fastq)

    # further filter these reads to find those matching the reverse replicate primer
    # Note that this only needs to be run in one sense since the reads have been rectified.
    rep_cutadapt_call = f"cutadapt -a {rep_primer}$ --discard-untrimmed -o {out_fastq} {for_temp_fastq}"
    run_sys_command(rep_cutadapt_call)    

    # If the spacer is not on the same strand as for_primer we need to reverse complement it in order to be directly analyzing A->G transisitons
    # and not C->T.
    if not spacer_for:
        print("Spacer and forward primer on opposite strands.")
        print("Reversing sequences.")
        out_fastq_temp = for_temp_fastq.with_name("out_fastq_temp.fastq")
        run_sys_command(f"seqtk seq -r {out_fastq} > {out_fastq_temp}")
        run_sys_command(f"mv {out_fastq_temp} {out_fastq}")
                        
    os.remove(for_temp_fastq)
    total_reads = get_fastq_read_count(out_fastq)
    print(f"Writing fastq file: {out_fastq} with {total_reads} entries.")
    return

def build_pileup(filtered_fastq, amplicon_fasta, thread_count = 10):
    """
    Make a pileup from a filtered fastq and the fasta of the amplicon.

    Usage:
    build_pileup(filtered_fastq, amplicon_fasta, thread_count = 10)

    Description:
    This function assigns stacked deletion variants to the fastq reads using the deletion matrix
    representation of the variants. It uses minimap2 to align the fastq reads to the WT reference.

    Inputs:
    filtered_fastq: str or path, path to filtered fastq file made with extract_rxn_fastq()
    amplicon_fasta: str or path, path to fasta file of amplicon sequence
    thread_count: int, number of threads to use a various stages that allow multi-threading
    
    Returns:
    df_pileup: dataframe with the pileup indexed by amplicon position
    """
    
    # make a bowtie2 compatible index for alignment
    bowtie_ind = Path(amplicon_fasta).with_name(f"{Path(amplicon_fasta).stem}_idx")
    bowtie_ref = f"bowtie2-build {amplicon_fasta} {bowtie_ind}"
    print(f"Generating bowtie2 reference for {amplicon_fasta}.")
    run_sys_command(bowtie_ref)
    print("")

    # align fastq reads to amplicon reference
    print(f"Aligning {filtered_fastq} to amplicon sequence.")
    sam_align = Path(filtered_fastq).with_name(f"{Path(filtered_fastq).stem}_align.sam")
    bowtie_align = f"bowtie2 -x {bowtie_ind} -U {filtered_fastq} -S {sam_align} -p {thread_count}"
    run_sys_command(bowtie_align)
    print("")

    # Get to a BAM file that can be used for the pileup
    print("Preparing binary form of alignment, ie. SAM to BAM.")
    bam_align = sam_align.with_name(f"{sam_align.stem}.bam")
    sam_sort_command = f"samtools view -bS {sam_align} -@ {thread_count} | samtools sort -o {bam_align} -@ {thread_count}"
    run_sys_command(sam_sort_command)
    run_sys_command(f"samtools index {bam_align}")
    print("")

    # Build the pileup with samtools
    print("Running pileup analysis and calculating edits.")
    pileup_file = bam_align.with_name(f"{Path(filtered_fastq).stem}.pileup")
    pileup_call = f"samtools mpileup -f {amplicon_fasta} {bam_align} > {pileup_file}"
    run_sys_command(pileup_call)
    print("")

    # parse the pileup file into a dataframe
    print("Reading pileup into dataframe.")
    columns = ["Reference", "Position", "RefBase", "Coverage", "Bases", "Qualities"]
    df_pileup = pd.read_csv(pileup_file, sep="\t", names=columns).set_index('Position')
    print("")

    # delete all the temporary files that were created
    print("Cleaning up temporary files.")
    
    for btind in glob(str(Path(amplicon_fasta).parent/'*.bt2')):
        os.remove(btind)
    
    os.remove(sam_align)
    os.remove(bam_align)
    os.remove(pileup_file)
    print("")
    
    print("Finished.")
    print("-----------------------------------------------------------------------")

    return df_pileup

def analyze_edits(df_pileup, positions, min_qscore=25):
    """
    Take the pileup dataframe and some positions of interest and summarizes the mutations at those positions 

    Usage:
    analyze_edits(df_pileup, positions, min_qscore=25)

    Inputs:
    df_pileup: dataframe, output from build_pileup()
    positions: list of ints, list of positions in the amplicon to analyze
    min_qscore: float, minimum q-score for base to be counted.
    
    Returns:
    edit_df: dataframe indexed by position with the read count, the WT base, and the percentage of occurance of each base
    """
    # initialize dictionary for edits
    edit_dict = {}
    # iterate over positions
    for position in positions:
        # initialize dictionary entry for position
        edit_dict[position] = {'read_count':0, 'WT_base':'', 'A': 0, 'G': 0, 'T': 0, 'C': 0, 'other':0}
        # fetch the pileup entry for position
        posrow = df_pileup.loc[position]
        # get the pileup strings into a manageable format
        editdf = parse_base_string(posrow['Bases'], posrow['Qualities'])
        # take only those reads that have the correct sense (which everything should given the processing)
        # and that meet the qscore threshold. Get the unique base counts.
        basecounts = editdf[(editdf['qscore'] > min_qscore) & (editdf['sense'] == '+')]['base'].value_counts()
        # get the total number of reads
        total_counts = basecounts.sum()
        # store the total number of reads in the entry
        edit_dict[position]['read_count']=total_counts
        # store the WT base
        edit_dict[position]['WT_base'] = posrow['RefBase']
        # iterate over all of the unique bases present
        for bc in basecounts.index:
            # if the base is '.' use the WT base since '.' stands for match to reference
            if bc == '.':
                key = posrow['RefBase']
            # if the base is among AGTC use it directly
            elif bc in 'AGTC':
                key = bc
            # if it is something else, classify it as other
            else:
                key = 'other'
            # put the percentage of reads matching the 'key' base into the dictionary entry
            edit_dict[position][key] = 100*(basecounts[bc] / total_counts)
    # turn the dictionary into a dataframe
    edit_df = pd.DataFrame.from_dict(edit_dict, orient='index')
    return edit_df

def process_tube_rxns(rxn_table0, primers_table0, tube, tube_fastq):
    """
    Iterates over the reactions in a given pooled tube and takes them to pileup dataframes. 

    Usage:
    process_tube_rxns(rxn_table0, primers_table0, tube, tube_fastq)

    Description:
    This function segregates all of the reads from a given pooled sequencing sample into their respective conditions according to
    the tabulated information in rxn_table0 and primers_table0. It aligns and performs a pileup analysis for all conditions separately
    and returns a dictionary containing all of the pileups. Each pileup can be analyzed for editing with analyze_edits()

    Inputs:
    rxn_table0: dataframe, from csv with reaction information, columns are: {'rxn', 'tube', 'for_primer', 'rev_primer'}. 'tube' matches the id of the reaction tube,
        'for_primer' is the id of the forward primer, 'rev_primer' is the id of the reverse (or replicate) primer
    primers_table0: dataframe, from csv with primer information, columns are: {'id' (as idx), 'seq', 'base_editor', 'site', 'sense', 'Amplicon', 'spacer_for'}
        'id' is the primer id used for lookup, 'seq' is the primer sequence, 'base_editor' is the variant encoded by the primer if forward, 'site' is the id of
        the amplicon site, 'sense' is the type of primer (either for for the base_editor barcode or rev for the replicate barcode), 'replicate' is the replicate id,
        'Amplicon' is the amplicon sequence in the sense of the forward primer, 'spacer_for' is a boolean for whether the spacer is on the same or opposite strand 
        to the forward primer.
    tube: int, tube number id
    tube_fastq: path to the consensus fastq (made by paired_to_consensus) that contains the reads for the tube
    
    Returns:
    pileup_dict: dict, a dictionary of all conditions (ie. variant_site_replicate) containing the pileup_df for the reads matching those conditions.
    """

    # initialize dictionary
    pileup_dict = {}
    # iterate over all of the rxn in the pooled tube
    for row in rxn_table0[rxn_table0['tube'] == tube].iterrows():
        # get the primer ids for the forward and reverse primers in the rxn
        # and use them to look up information in the primer_table0 for sequences, sites, variants, replicates, etc.
        for_primer_id = row[1]['for_primer']
        rev_primer_id = row[1]['rev_primer']
        for_primer = primers_table0.loc[for_primer_id,'seq']
        rev_primer = primers_table0.loc[rev_primer_id,'seq']
        base_editor = primers_table0.loc[for_primer_id,'base_editor']
        amp_site_for = primers_table0.loc[for_primer_id,'site']
        amp_site_rev = primers_table0.loc[rev_primer_id,'site']
        # Check that the amplicon sites agree between the forward and reverse primers
        if amp_site_for != amp_site_rev:
            raise Exception(ValueError, "Amplicon sites don't match!")
        # generate the name for the reaction to specify the base_editor variant, the amplicon site, and the replicate number.
        rxn_name = f"{base_editor}_amp{amp_site_for}_rep{primers_table0.loc[rev_primer_id,'replicate']}"
        print(rxn_name)
        extract_fastq = Path(tube_fastq).with_name(f"{rxn_name}.fastq")

        spacer_forward = primers_table0.loc[for_primer_id,'spacer_for']
        # run the extraction script to filter out reads matching the rxn conditions.
        extract_rxn_fastq(tube_fastq, for_primer, rev_primer, extract_fastq, spacer_for=spacer_forward)
        
        # write amplicon fasta but need to check if the spacer is opposite the forward primer, in which case the amplicon must be reverse complemented.
        amplicon_seq = Seq.Seq(primers_table0.loc[for_primer_id,'Amplicon'])
        
        if not spacer_forward:
            amplicon_seq = amplicon_seq.reverse_complement()
                               
        # build the amplicon sequence record
        amplicon_record = SeqRecord.SeqRecord(seq = amplicon_seq,
                                              id = f"amp{amp_site_for}",
                                              name = f"amp{amp_site_for}",
                                              description = '')
        amplicon_fasta = extract_fastq.with_name(f"amp{amp_site_for}.fasta")
        # write the amplicon fasta to file
        SeqIO.write(amplicon_record, amplicon_fasta, 'fasta')
        # run build_pileup on the extracted fastq and add its pileup dataframe entry to dictionary under the rxn name key
        pileup_dict[rxn_name] = build_pileup(extract_fastq, amplicon_fasta, thread_count = 10)
        
    return pileup_dict

def evaluate_adenines(pileup_dict):
    """
    Take the pileup dictionary and evaluate the edits at all adenines. 

    Usage:
    evaluate_adenines(pileup_dict)

    Inputs:
    pileup_dict: dict, a dictionary of all conditions (ie. variant_site_replicate) containing the pileup_df for the reads matching those conditions.
    
    Returns:
    aedit_dict: dict, a dictionary of all conditions (ie. variant_site_replicate) containing the edit_df for the reads matching those conditions.    
    """
    aedit_dict = {}
    for key in list(pileup_dict.keys()):
        pileup = pileup_dict[key]
        pileupA = pileup[pileup['RefBase'] == 'A']
        pileupA_df = analyze_edits(pileupA, pileupA.index)
        print(key)
        print(f"Maximum editing: {pileupA_df['G'].max()}%")
        print('')
        aedit_dict[key] = pileupA_df
    return aedit_dict