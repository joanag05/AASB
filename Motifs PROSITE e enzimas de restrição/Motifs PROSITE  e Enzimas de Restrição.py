class Rest_Enzy:
    """Classe para análise de enzimas de restrição
    """

    def __init__(self, enzyme_seq, seq):
        """Construtor da classe

        Args:
            enzyme_seq (str): Sequência de uma enzima de restrição com local de corte
            seq (str): Sequência de DNA
        """
        self.enzyme_seq = enzyme_seq
        self.seq = seq

    
    def __str__(self):
        """Representação das sequências
        """
        return f'The sequence we want to study is {self.enzyme_seq} and the restriction enzyme is {self.seq}'


    def __repr__(self):
        """Representação das sequências 
        """
        return f'Rest_Enzy(\'{self.enzyme_seq}\', {self.seq})'


    def Regex(self): 
        """Retira o local de corte da sequência da enzima de restrição e possibilita as ambiguidades 

        Returns:
            str: Expressão regular 
        """
        dic = {"A":"A", "C":"C", "G":"G", "T":"T", "R":"[GA]", "Y":"[CT]", "M":"[AC]", "K":"[GT]",
         "S":"[GC]", "W": "[AT]", "B":"[CGT]", "D":"[AGT]", "H":"[ACT]", "V":"[ACG]", "N":"[ACGT]"}                
        corte = self.enzyme_seq.replace("^","")
        regex = ""
        for x in corte:
            regex += dic[x]
        return regex   


    def positions(self):  
        from re import finditer
        """Dada uma enzima de restrição e uma sequência de DNA indica, se possível, as posições do corte

        Returns:
            list: Lista com as posições de corte. Se não existir corte retorna lista vazia
        """
        cutpos = self.enzyme_seq.find("^")
        regex = Rest_Enzy.Regex(self)
        matches = finditer(regex, self.seq)
        locs = []
        for x in matches:
            locs.append(x.start() + cutpos)
        return locs   
 

    def subsequences (self):
        """A partir das posições de corte permite obter as subsequências geradas pelo corte da enzima

        Returns:
            list: Lista com as subsquências geradas pelo corte da enzima. No caso de não haver corte retorna 
                uma lista com um elemento único que representa a sequência original
        """
        res = []
        positions = Rest_Enzy.positions(self)
        positions.insert(0,0)
        positions.append(len(self.seq))
        for i in range(len(positions)-1):
            res.append(self.seq[positions[i]:positions[i+1]])
        return res   


    def positions_rc(self):
        """Origina a sequência complementar reversa e procura por locais de corte
        Returns:
            list: Lista com as posições de corte (no sentido 5' -> 3'). Se não existir corte retorna lista vazia
        """
        from re import finditer            
        temp = self.seq.replace("A","t").replace("T","a").replace("C","g").replace("G","c")
        rev_comp = temp[::-1].upper()
        cutpos = self.enzyme_seq.find("^")
        regex = Rest_Enzy.Regex(self)
        matches = finditer(regex, rev_comp)
        locs_rc = []
        for x in matches:
            locs_rc.append(x.start() + cutpos)
        return locs_rc     


    def subsequences_rc(self):
        """A partir das posições de corte permite obter as subsequências geradas pelo corte da enzima na 
        cadeia complementar reversa

        Returns:
            list: Lista com as subsquências geradas pelo corte da enzima. No caso de não haver corte retorna 
                uma lista com um elemento único que representa a sequência original
        """
        res = []
        positions = Rest_Enzy.positions_rc(self)
        positions.insert(0,0)
        positions.append(len(self.seq))
        temp = self.seq.replace("A","t").replace("T","a").replace("C","g").replace("G","c")
        rev_comp = temp[::-1].upper()
        for i in range(len(positions)-1):
            res.append(rev_comp[positions[i]:positions[i+1]])
        return res

import unittest

class Test_Rest_Enzy(unittest.TestCase):
    """Classe para a testagem da class Rest_Enzy
    """

    def test_Regex(self):
        """Teste ao método Regex
        """
        self.assertEqual(Rest_Enzy("G^AMTV", "GTAGAAGATTCTGACTGATCGATTC").Regex(), "GA[AC]T[ACG]")
        self.assertEqual(Rest_Enzy("G^AATC", "GTAGAAGATTCTGACTGATCGATTC").Regex(), "GAATC")
        self.assertEqual(Rest_Enzy("GACT^A", "GTAGAAGATTCTGACTGATCGATTC").Regex(), "GACTA") 
        self.assertEqual(Rest_Enzy("A^TRCGMA", "GTAGATACGCAGATTCTGATGCGCACT").Regex(), "AT[GA]CG[AC]A")


    def test_positions(self):
        """Teste ao método positions
        """
        self.assertEqual(Rest_Enzy("G^AMTV", "GTAGAAGATTCTGACTGATCGATTC").positions(), [13])
        self.assertEqual(Rest_Enzy("G^AATC", "GTAGAAGATTCTGACTGATCGATTC").positions(), [])
        self.assertEqual(Rest_Enzy("GACT^A", "GTAGAAGATTCTGACTGATCGATTC").positions(), []) 
        self.assertEqual(Rest_Enzy("A^TRCGMA", "GTAGATACGCAGATTCTGATGCGCACT").positions(), [5, 19])


    def test_subsequences(self):
        """Teste ao método subsequences
        """
        self.assertEqual(Rest_Enzy("G^AMTV", "GTAGAAGATTCTGACTGATCGATTC").subsequences(), ['GTAGAAGATTCTG', 'ACTGATCGATTC'])
        self.assertEqual(Rest_Enzy("G^AATC", "GTAGAAGATTCTGACTGATCGATTC").subsequences(), ['GTAGAAGATTCTGACTGATCGATTC'])
        self.assertEqual(Rest_Enzy("GACT^A", "GTAGAAGATTCTGACTGATCGATTC").subsequences(), ['GTAGAAGATTCTGACTGATCGATTC']) 
        self.assertEqual(Rest_Enzy("A^TRCGMA", "GTAGATACGCAGATTCTGATGCGCACT").subsequences(), ['GTAGA', 'TACGCAGATTCTGA', 'TGCGCACT'])    

    
    def test_positions_rc(self):
        """Teste ao método positions_rc
        """
        self.assertEqual(Rest_Enzy("G^AMTV", "GTAGAAGATTCTGACTGATCGATTC").positions_rc(), [1, 15])
        self.assertEqual(Rest_Enzy("G^AATC", "GTAGAAGATTCTGACTGATCGATTC").positions_rc(), [1, 15])
        self.assertEqual(Rest_Enzy("TC^AGTC", "GTAGAAGATTCTGACTGATCGATTC").positions_rc(), [9]) 
        self.assertEqual(Rest_Enzy("A^TRCGMA", "GTAGATACGCAGATTCTGATGCGCACT").positions_rc(), [])

    
    def test_subsequences_rc(self):
        """Teste ao método subsequences_rc
        """
        self.assertEqual(Rest_Enzy("G^AMTV", "GTAGAAGATTCTGACTGATCGATTC").subsequences_rc(), ['G', 'AATCGATCAGTCAG', 'AATCTTCTAC'])
        self.assertEqual(Rest_Enzy("G^AATC", "GTAGAAGATTCTGACTGATCGATTC").subsequences_rc(), ['G', 'AATCGATCAGTCAG', 'AATCTTCTAC'])
        self.assertEqual(Rest_Enzy("GACT^A", "GTAGAAGATTCTGACTGATCGATTC").subsequences_rc(), ['GAATCGATCAGTCAGAATCTTCTAC']) 
        self.assertEqual(Rest_Enzy("A^TRCGMA", "GTAGATACGCAGATTCTGATGCGCACT").subsequences_rc(), ['AGTGCGCATCAGAATCTGCGTATCTAC'])


unittest.main(argv=[''], exit=False)