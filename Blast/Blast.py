
import mysql.connector

db = mysql.connector.connect(host="127.0.0.1",   
                     user="root",         
                     passwd="Cardoso10@",  
                     db="trabalho",
                     auth_plugin = "mysql_native_password"
                     ) 


def obter_seq(host, user, passwd, db, auth_plugin):
    
    cur = db.cursor()
    cur.execute("""
        SELECT Seq from Sequence
        join Gene_Bank on Sequence.ID = Gene_Bank.Keywords
        where organism="Candida albicans"
        """)
    sequence = cur.fetchall()
    
    return sequence[0]


def obter_seqs_BD(host, user, passwd, db, auth_plugin):

    cur = db.cursor()
    cur.execute("""
        SELECT Seq from Sequence
        join Gene_Bank on Sequence.ID = Gene_Bank.Keywords
        """)
    sequences = cur.fetchall()
    return sequences 



class Blast:
    """Implementação de uma versão simplidicada do Blast 
    """

    def __init__(self, query, seq, w=3):
        """Construtor da classe

        Args:
            query (str): _description_
            seq (str): _description_
            w (int, optional): _description_. Defaults to 3.
        """
        self.query = query
        self.seq = seq
        self.w = w


    def query_map(query, w):
        """Obtém todas as palavras de tamanho w da query

        Args:
            query (str): query
            w (int): window size

        Returns:
            dict: Dicinário com chaves como subsquências de tamanho w e valores como offsets
        """
        tam = len(query)
        res = {}
        for chave, offset in [(query[p: p + w], p) for p in range(0, tam - w + 1)]:
            if chave not in res: res[chave] = []
            res[chave].append(offset)
        return res


    def get_all_offsets(s1, s2):
        """Determinação dos offsets
        """
        w = len(s1)
        res = []
        for p in range(0, len(s2) - w + 1):
            if s2[p : p + w] == s1:
                res.append(p)
        return res


    def hits(qm, seq):
        """ Recebe o dicionário da função query_map e uma sequência da BD
        Args:
            qm (dict): dicinário previamente gerado com chaves como subsquências de tamanho w e valores como offsets
            seq (str): sequência da BD

        Returns:
            list: Lista de hits em que cada elemento é um tuplo com os índices
        """
        res = []
        for chave, offsets in qm.items():
            for o_query in offsets:
                for o_seq in Blast.get_all_offsets(chave, seq):
                    res.append((o_query, o_seq))
        return res


    def extend_hit_dir(query, seq, o1, o2, direction):
        """Função auxiliar que permite extender o hit em ambas as direções
        """
        matches = 0
        count = 0
        while o1 >= 0 and o2 >= 0 and o1 < len(query) and o2 < len(seq):
            matches += 1 if query[o1] == seq[o2] else 0
            count += 1
            if 2 * matches < count:
                return o1, o2, matches, count - 1
            o1 += direction
            o2 += direction
        return o1 - direction, o2 - direction, matches, count


    def extend_hit(query, seq, hit, w):
        """recebe a query, a sequência da BD, o hit e o valor de w e o estende um hit em cada
        direção se o nº de matches for de pelo menos metade do tamanho da extensão;

        Args:
            query (str): query
            seq (str): sequência da BD
            hit (tuple): Hit
            w (int): window size

        Returns:
            tuple: tuplo com o índice do início do hit
        """
        o1, o2 = hit
        left  = Blast.extend_hit_dir(query, seq, o1 - 1, o2 - 1, -1)
        right = Blast.extend_hit_dir(query, seq, o1 + w, o2 + w, +1)

        O1, O2, ML, SL = left
        _,   _, MR, SR = right

        return O1, O2, w + SL + SR, ML + w + MR
        

    def best_hit(query, seq, w):
        """recebe a query, a sequência da BD e o valor de w e calcula o melhor hit 

        Args:
            query (str): query
            seq (str): sequência da BD
            w (int): window size

        Returns:
            tuple: tuplo correspondente ao melhor hit
        """
        hits = Blast.hits(Blast.query_map(query, w), seq)
        bestScore = -1.0
        best = ()
        for h in hits:
            ext = Blast.extend_hit(query, seq, h, w)
            score = ext[3]
            if score > bestScore or (score== bestScore and ext[2] < best
        [2]):
                bestScore = score
                best = ext
        return best

import unittest


class teste_Blast(unittest.TestCase):
    """Classe para a testagem da class Blast
    """

    def test_query_map(self):
        """Teste ao método query_map
        """
        self.assertEqual(Blast.query_map("AATATAT", w = 3), {'AAT': [0], 'ATA': [1, 3], 'TAT': [2, 4]})

    
    def test_get_all_offsets(self):
        """Teste ao método get_all_offsets
        """
        self.assertEqual(Blast.get_all_offsets("AATATAT","AATATGTTATATAATAATATTT"), [])


    def test_hits(self):
        """Teste ao método hits 
        """
        self.assertEqual(Blast.hits(Blast.query_map("AATATAT", 3), "AATATGTTATATAATAATATTT"),[(0, 0), (0, 12), (0, 15), (1, 1), (1, 8), (1, 10),
                                                                    (1, 13), (1, 16), (3, 1), (3, 8), (3, 10),
                                                                    (3, 13), (3, 16), (2, 2), (2, 7), (2, 9), (2, 17),
                                                                    (4, 2), (4, 7), (4, 9),(4, 17)])
        
    
    def test_extend_hit_dir(self):
        """Teste ao método extend_hit_dir 
        """
        self.assertEqual(Blast.extend_hit_dir("AATATAT", "AATATGTTATATAATAATATTT", 1, 16, -1), (0, 15, 2, 2))

    
    def test_extend_hit(self):
        """Teste ao método extend_hit
        """
        self.assertEqual(Blast.extend_hit("AATATAT", "AATATGTTATATAATAATATTT", (1,16), 3), (0, 15, 7, 6))


    def test_best_hit(self):
        """Teste ao método best_hit
        """
        self.assertEqual(Blast.best_hit("AATATAT", "AATATGTTATATAATAATATTT", 3), (0, 0, 7, 6))

        
unittest.main(argv=[''], exit=False)