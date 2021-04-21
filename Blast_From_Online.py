from PrimerBlast_Tools import Primer_Blast_Fuction as pb

gene_id = {
    'db': 'nucest',
    'id': 'NM_120759.3',
    'rettype': 'gb',
    'retmode': 'text'
}
email = 'sample@******.com'
gene_name = 'AtLEA4-5'
out_address = '/Users/******/Desktop/******/'         # Needs a '/' in the end
db_address = '/Users/******/Desktop/ncbi-blast-2.10.1+/bin/******.fa'
print('Gene information recieved.')
evalue = 1
identity = 10

record, target_seq = pb.get_target(email, gene_id, gene_name)
print('\n', record, '\n')
print(len(target_seq))

count,name=pb.blastn(query_address=gene_name, db_address=db_address, out_address1=out_address + gene_name + '_Blast.xml', evalue=evalue, identity=identity, task='megablast', dust='yes')
