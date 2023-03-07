from Bio.Blast import NCBIXML

# Parse the results from BlastP
result_handle = open("blastp.xml")
blast_records = NCBIXML.parse(result_handle)

# Iterate over each Blast record
for blast_record in blast_records:
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            # Get the alignment sequences
            query_seq = hsp.query
            subject_seq = hsp.sbjct
            # Initialize variables to keep track of the mismatches
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
            # Print the results for this HSP
            print(f"Query: {blast_record.query}")
            #print("HSP name:", blast_record.query, alignment.hit_def, alignment.title)
            print("HSP name:",alignment.title)
            #print(f"Alignment length: {hsp.align_length}")
            #print(f"Mismatches: {mismatches}")
            #print("Mismatch positions:")
            for pos in mismatch_positions:
                #print(f"Position {pos[0]}: Query AA = {pos[1]}, Ref AA = {pos[2]}")
                print(pos[1]+str(pos[0])+ pos[2])
            print()

