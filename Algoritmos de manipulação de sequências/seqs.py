class SEQS:
    """class para análise de sequências biológicas
    """
    def __init__(self, seq):
        """Construtor da class

        Args:
            seq (str): Sequência de DNA
        """
        self.seq = seq

    def __str__(self):
        return f'A sequência a analisar é: {self.seq}'
    
    def __repr__(self):
        return f'SEQS(\"{self.seq}\")'
    
    def __getitem__(self,n):
        return self.seq[n]

    def __len__(self):
        return len(self.seq)


    def biotype(self):
        """ Função recebe uma string

        Retorna:
            booleano: DNA, RNA, AMINO ou ERRO
        """
        if len([x for x in self.seq if x not in  "ACGTactg"]) == 0:
            return ('DNA')
        elif len([x for x in self.seq if x not in  "ACGU"]) == 0:
            return ('RNA')
        elif len([x for x in self.seq if x not in  "ABCDEFGHIJKLMNOPQRSTUVWXYZ"]) == 0:
            return('AMINO')
        else:
            return("ERRO")


    def count_nucleo(self):
        """Função que recebe uma sequência de DNA
        
        Retorna:
            list: lista com o número de A, C, G e T
        """
        A=self.seq.count('A')
        C=self.seq.count('C')
        G=self.seq.count('G')
        T=self.seq.count('T')
        return [f"A ---> {A} C ---> {C} G ---> {G} T ---> {T}"]


    def freq (self):
        """Função que calcula a frequência de cada símbolo na sequência.

        Retorna:
            Um dicionário com a frequência de cada núcleótido
        """
        dic = {}
        for s in self.seq.upper():
            if s in dic: dic[s] += 1
            else : dic[s] = 1
        return dic

    def gc_content(self):
        """Função que calcula percentagem de nucleotídeos G e C numa sequência de DNA

        Retorna:
            float: A percentagem de nucleotídeos G e C numa sequência de DNA
        """
        gc_count = 0
        for o in self.seq:
            if o in "GCgc": gc_count += 1
        return gc_count / len (self.seq)

    

    def transcription(self):
        """ Função que calcula o RNA correspondente á transcrição da sequência de DNA fornecida.

        Retorna:
            str: A sequência de RNA correspondente á transcrição da sequência de DNA fornecida. 
        """
        if(SEQS(self.seq).biotype() !=  'DNA') : return 'NOT DNA'
        return self.seq.upper().replace('T','U') 
        

    
    def reverse_complement(self):
        if(SEQS(self.seq).biotype() !=  'DNA') : return 'NOT DNA'
        """Função que recebe uma sequência de DNA 


        Returna:
            str:O complemento reverso da sequência de DNA 
        """
        self.dna_inv = self.seq[:: -1].upper()
        self.dna_inv = self.dna_inv.replace('A', 't').replace('T', 'a').replace('G', 'c').replace('C', 'g')
        return self.dna_inv 

    def get_codons(self):
        if(SEQS(self.seq).biotype() !=  'DNA') : return 'NOT DNA'
        """Função que recebe uma sequência de DNA 

        Retorna:
            list: lista de codões
        """
        DNA =self.seq.upper()
        codons = []
        for i in range(0, len(DNA), 3):
            codons.append(DNA[i : i + 3])
        if len(codons[-1]) < 3:
            codons.pop()
        return codons

    def codon_to_amino(self):
        """função que recebe a lista de codões

        Retorna:
            str: sequência de aminoácidos
        """
        gencode = {
            'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
            'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
            'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
            'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
            'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
            'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
            'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
            'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
            'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
            'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
            'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
            'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
            'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
            'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
            'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
            'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'}
        
        if(SEQS(self.seq).biotype() !=  'DNA') : return 'NOT DNA'
        self.seq = self.seq.upper()
        return ''.join([gencode.get(self[p: p + 3],'')for p in range(0,len(self),3)])    
    
    def get_prots (self):
        """Função que recebe uma sequência de aminoácidos

        Retorna:
            list: lista de possíveis proteínas
        """
        if(SEQS(self.seq).biotype() !=  'DNA') : return 'NOT DNA'
        self.seq = self.seq.upper()
        proteins=[]   
        stop_codons_1 = SEQS(self.seq).codon_to_amino()
        stop_codons=stop_codons_1.split("_")
        for x in stop_codons:     
            for idx, amino_s in enumerate(x):    
                if amino_s == "M":                       
                    proteins.append(x[idx:]) 
        return proteins 
    
    def get_orfs (self):
        if(SEQS(self.seq).biotype() !=  'DNA') : return 'NOT DNA'
        """Função que recebe uma sequência de DNA

        Returns:
            list: lista com as seis ORFs
        """  
        DNA =self.seq.upper()
        orfs =[]
        orfs.append(DNA[0:])     
        orfs.append(DNA[1:])   
        orfs.append(DNA[2:])
        seq_inv = SEQS(self.seq).reverse_complement()
        orfs.append(seq_inv[0:])
        orfs.append(seq_inv[1:])  
        orfs.append(seq_inv[2:])
        return orfs
    
    def get_all_prots (self):
        """Função que recebe as 6 orfs

        Returns:
            list: lista com todas as proteínas
        """
                 
        orfs1 = SEQS(self.seq).get_orfs()
        codons1 = [] 
        aminos1 = []
        aminofinal=[]
        prots = []
        

        for idx, x in enumerate(orfs1):
            codons1.append(SEQS(x).get_codons())
    
        for y in (codons1):
            for b in y:
                aminos1.append(SEQS(b).codon_to_amino())
            aminofinal.append(''.join(aminos1))
            aminos1=[]
        
        for idx , k in enumerate(aminofinal):
                x = aminofinal[idx].split("_")
                try:
                    while True:
                        x.remove("")
                except ValueError:
                    pass
                
                for c in x:
                    for idx1, j in enumerate(c):  
                        if j == "M":                       
                            prots.append(c[idx1:])
        return prots




import unittest

class Test_SEQS(unittest.TestCase):
    """Classe que permite testagem do código do alinhamento NW
    """
    def test_transcription(self):
        """Testar a função transcription
        """
        self.assertEqual(SEQS("ATGATACGATGGCGTATAA").transcription(),'AUGAUACGAUGGCGUAUAA')
        self.assertEqual(SEQS("actgcata").transcription(),'ACUGCAUA')
        self.assertEqual(SEQS('JACGATAGCHKSS').transcription(),'NOT DNA')
        self.assertEqual(SEQS('MIRWRI').transcription(),'NOT DNA')

    def test_reverse_complement(self):
        """Testar a função reverse_complement
        """
        self.assertEqual(SEQS("ATGATACGATGGCGTATAA").reverse_complement(),'ttatacgccatcgtatcat')
        self.assertEqual(SEQS("actgcata").reverse_complement(),'tatgcagt')
        self.assertEqual(SEQS('JACGATAGCHKSS').reverse_complement(),'NOT DNA')
        self.assertEqual(SEQS('MIRWRI').reverse_complement(),'NOT DNA')

    def test_get_codons(self):
        """Testar a função get_codons
        """
        self.assertEqual(SEQS("ATGATACGATGGCGTATAA").get_codons(),['ATG', 'ATA', 'CGA', 'TGG', 'CGT', 'ATA'])
        self.assertEqual(SEQS("actgcata").get_codons(),['ACT', 'GCA'])
        self.assertEqual(SEQS('JACGATAGCHKSS').get_codons(),'NOT DNA')
        self.assertEqual(SEQS('MIRWRI').get_codons(),'NOT DNA')

    def test_codon_to_amino(self):
        self.assertEqual(SEQS("ATGATACGATGGCGTATAA").codon_to_amino(),'MIRWRI')
        self.assertEqual(SEQS("actgcata").codon_to_amino(),'TA')
        self.assertEqual(SEQS('JACGATAGCHKSS').codon_to_amino(),'NOT DNA')
        self.assertEqual(SEQS('MIRWRI').codon_to_amino(),'NOT DNA')

    def test_get_prots(self):
        self.assertEqual(SEQS("ATGATACGATGGCGTATAA").get_prots(),['MIRWRI'])
        self.assertEqual(SEQS('agtagaatgattg').get_prots(),['MI'])
        self.assertEqual(SEQS('JACGATAGCHKSS').get_prots(),'NOT DNA')
        self.assertEqual(SEQS('MIRWRI').get_prots(),'NOT DNA')

    def test_get_orfs(self):
        self.assertEqual(SEQS("ATGATACGATGGCGTATAA").get_orfs(),['ATGATACGATGGCGTATAA','TGATACGATGGCGTATAA','GATACGATGGCGTATAA','ttatacgccatcgtatcat','tatacgccatcgtatcat','atacgccatcgtatcat'])
        self.assertEqual(SEQS('agtagaatgattg').get_orfs(),)

       

        







        

    

   
if __name__ == '__main__':
    unittest.main()