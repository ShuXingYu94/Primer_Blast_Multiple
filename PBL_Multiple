from PrimerBlast_Tools import Primer_Blast_Fuction as pb

work_direction = '/Users/******/Desktop/******.fa'    # This works fine with txt or fasta file
out_address = '/Users/******/Desktop/******/'         # Needs a '/' in the end
db_address = '/Users/******/Desktop/ncbi-blast-2.10.1+/bin/******.fa'
evalue = 1
identity = 18
FASTA = open(work_direction)
global_args = {
    'PRIMER_OPT_SIZE': 20,
    'PRIMER_PICK_INTERNAL_OLIGO': 0,
    'PRIMER_INTERNAL_MAX_SELF_END': 8,
    'PRIMER_MIN_SIZE': 18,
    'PRIMER_MAX_SIZE': 25,
    'PRIMER_OPT_TM': 60.0,
    'PRIMER_MIN_TM': 50.0,
    'PRIMER_MAX_TM': 68.0,
    'PRIMER_MIN_GC': 20.0,
    'PRIMER_MAX_GC': 80.0,
    'PRIMER_MAX_POLY_X': 100,
    'PRIMER_INTERNAL_MAX_POLY_X': 100,
    'PRIMER_SALT_MONOVALENT': 50.0,
    'PRIMER_DNA_CONC': 50.0,
    'PRIMER_MAX_NS_ACCEPTED': 0,
    'PRIMER_MAX_SELF_ANY': 12,
    'PRIMER_MAX_SELF_END': 8,
    'PRIMER_PAIR_MAX_COMPL_ANY': 12,
    'PRIMER_PAIR_MAX_COMPL_END': 8,
    'PRIMER_PRODUCT_SIZE_RANGE': [[100, 500], [500, 1000]],
    'PRIMER_NUM_RETURN': 10
}

print('Parameters and Sequence(s) Recieved.')

seqs, names = pb.readFASTA(FASTA)

for a in range(0,len(names)):
    raw_seq=seqs[a]
    gene_name=names[a]
    print('\n'+'-'*10+'Start designing primers for: ' + gene_name+ '.'+'-'*10)

    raw_seq = raw_seq.upper()
    Nuc = list('ATGC')
    amb_Nuc = list('NRYKMSWBDHV')
    target_seq = ''
    for nuc in raw_seq:
        if nuc in amb_Nuc:
            print('Error: Found Ambiguous DNA base.')
            exit(1)
        if nuc == 'U':
            target_seq += 'T'
        if nuc in Nuc:
            target_seq += nuc

    seq_args = {  # Set parameters: name,seq, include region etc.
        'SEQUENCE_ID': gene_name,
        'SEQUENCE_TEMPLATE': target_seq,
        'SEQUENCE_INCLUDED_REGION': [0, len(target_seq)],
        'SEQUENCE_TARGET': [],
        'SEQUENCE_EXCLUDED_REGION': []
    }
    primer_df = pb.design_primer(seq_args, global_args, out_address + gene_name + '_Original_Primer.txt')
    # To get all primer pairs in the form of dataframe
    primer_list = pb.extract_pairs(out_address + gene_name + '_Original_Primer.txt')
    # Turn dataframe into list of primer pairs
    print('\nStart confirming primer pairs...')

    list_result = []
    count = 0
    list_r = []
    list_f = []

    for pair in primer_list:
        print('\nChecking Primer pair ', str(count + 1), ':')
        forward = pair[0]
        tmp_f = open(out_address + "tmp_f.fasta", "w")
        tmp_f.write("%s" % forward)
        tmp_f.close()
        reverse = pair[1]
        tmp_r = open(out_address + "tmp_r.fasta", "w")
        tmp_r.write("%s" % reverse)
        tmp_r.close()
        print('_Left Primer:')
        f, f_n = pb.blastn(out_address + "tmp_f.fasta", db_address=db_address,
                           out_address1=out_address + gene_name + '_PBL.xml',
                           evalue=evalue, identity=len(forward) - 2, task='blastn-short', dust='no')
        print('_Right Primer:')
        r, r_n = pb.blastn(out_address + "tmp_r.fasta", db_address=db_address,
                           out_address1=out_address + gene_name + '_PBL.xml',
                           evalue=evalue, identity=len(reverse) - 2, task='blastn-short', dust='no')
        list_f.append(f)
        list_r.append(r)

        if r == 1 and f == 1:
            list_result.append(count)
            # print('Primer Pair '+str(count)+': Forward/Reverse primer found '+str(f)+'/'+str(r)+' products on intended
            # targets,respectively.')
        count += 1

    primer_df.loc['LEFT_PRIMER_PRODUCTS'] = list_f
    primer_df.loc['RIGHT_PRIMER_PRODUCTS'] = list_r
    pd_result = primer_df.iloc[:, list_result]
    if pd_result.empty:
        print('\nNo specified primer pair found. Try change parameters or increase primer numbers to return.\n')
    else:
        print('\nResult of useful primer(s):\n', pd_result)
        pd_result.to_csv(out_address + gene_name + '_PrimerBlast_Result.csv')
        print('\nResult saved.')
    for path in ['_Original_Primer.txt','_PBL.xml']:
        if pb.os.path.exists(out_address + gene_name + path):
            pb.os.remove(out_address + gene_name + path)

for path in ["tmp_f.fasta", "tmp_r.fasta"]:
    if pb.os.path.exists(out_address + path):
        pb.os.remove(out_address + path)
