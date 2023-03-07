from Bio.Blast import NCBIXML

# Parse the results from BlastP
result_handle = open("blastp.xml")
blast_records = NCBIXML.parse(result_handle)

# Initialize variables to keep track of the mismatches
mismatches = 0
mismatch_positions = []

# Iterate over each Blast record
for blast_record in blast_records:
    # Initialize variables for this query
    best_hsp = None
    best_evalue = float('inf')
    # Iterate over each alignment for this query
    for alignment in blast_record.alignments:
        # Iterate over each HSP for this alignment
        for hsp in alignment.hsps:
            # Select the HSP with the lowest E-value for this query
            if hsp.expect < best_evalue:
                best_hsp = hsp
                best_evalue = hsp.expect

# Get the alignment sequences for the best HSP for this query
    query_seq = best_hsp.query
    subject_seq = best_hsp.sbjct
    # Reset the variables for the mismatches for this query
    mismatches = 0
    mismatch_positions = []
    # Iterate over each character in the alignment
    for i in range(len(query_seq)):
        # Check if the characters don't match
        if query_seq[i] != subject_seq[i]:
            mismatches += 1
            # Get the amino acid in the query and the reference in each mismatch position
            mismatch_query_aa = query_seq[i]
            mismatch_ref_aa = subject_seq[i]
            # Add the position and the amino acids to the mismatch_positions list
            mismatch_positions.append((i+1, mismatch_query_aa, mismatch_ref_aa))
    # Print the results for the best HSP for this query
    print(f"Query: {blast_record.query}")
    print("HSP name:",alignment.title)
    #print(f"Alignment length: {best_hsp.align_length}")
    #print(f"Mismatches: {mismatches}")
    #print("Mismatch positions:")
    for pos in mismatch_positions:
        #print(f"Position {pos[0]}: Query AA = {pos[1]}, Ref AA = {pos[2]}")
        print(pos[1]+str(pos[0])+pos[2])
    print()

