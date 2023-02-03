import io

class Mat:

    def __init__(self, rows, cols):
        self.mat = [[0 for c in range(cols)] for r in range(rows)]

    def numRows (self): return len(self.mat)

    def numCols (self): return len(self.mat[0])

    def __str__(self):  
        return '\n'.join(' '.join(str(val) for val in row) for row in self.mat)

    def __repr__(self):
        return str(self)

    def __getitem__ (self, n):
        return self.mat[n]

class Motifs:
    """Classe para análise de motifs. Permite a visualização das matrizes de frequências absolutas, pwm e pssm.
    Assim como o cálcula da sequência mais provável e a probabilidade de gerar uma determinada sequência
    """

    def __init__(self, seqs): 
        """Construtor da classe

        Args:
            seqs (str): String com várias sequências de DNA ou proteína
        """
        self.seqs = seqs
        self.alfabeto =  ""
        temp = set("".join(self.seqs))
        temp = sorted(temp)
        for x in temp:
            self.alfabeto += x


    def create_matrix(self):
        """ Recebe uma string com motifs

        Returns:
            list: A matriz das frequências absolutas dos respetivos caractéres introduzidos
        """
        import io
        n_colunas = len(self.seqs[0])
        n_linhas = len(self.alfabeto)
        M = []

        for l in range(0, n_linhas):
            M.append([0] * n_colunas)
        
        for x in self.seqs:                                  
            for i in range(len(x)):
                linha = self.alfabeto.index(x[i])
                M[linha][i] +=1

        Alfabeto = []
        for z in self.alfabeto:
            Alfabeto.append(z)
        zipped = list(zip(Alfabeto, M))
        for x, y in zipped:
            print(x, *y, sep= "  ")
            

    def PWM(self, pseudo = 1):
        """Recebe o valor de pseudo-contagens

        Args:
            pseudo (int, optional): Pseudo-contagem. Defaults to 1.

        Returns:
            list: Matriz PWM
        """
        import io
        n_colunas = len(self.seqs[0])
        n_linhas = len(self.alfabeto)
        M = []

        for l in range(0, n_linhas):
            M.append([0] * n_colunas)
        
        for x in self.seqs:                                  
            for i in range(len(x)):
                linha = self.alfabeto.index(x[i])
                M[linha][i] +=1
        pwm = M
        for l in range(len(pwm)):
            for c in range(len(pwm[l])):
                pwm[l][c] = round((pwm[l][c] + pseudo) / (len(self.seqs) + (pseudo * len(self.alfabeto))),2)    

        Alfabeto = []
        for z in self.alfabeto:
            Alfabeto.append(z)
        zipped = list(zip(Alfabeto, pwm))
        for x, y in zipped:
            print(x, *y, sep= "  ")


    def PSSM(self, pseudo = 0.5):
        """Recebe o valor de pseudo-contagens

        Args:
            pseudo (float, optional): _description_. Defaults to 0.5.

        Returns:
            list: Matriz PSSM
        """
        import math
        import io
        n_colunas = len(self.seqs[0])
        n_linhas = len(self.alfabeto)
        M = []

        for l in range(0, n_linhas):
            M.append([0] * n_colunas)
        
        for x in self.seqs:                                  
            for i in range(len(x)):
                linha = self.alfabeto.index(x[i])
                M[linha][i] +=1
        
        pssm = M
        for l in range(len(pssm)):
            for c in range(len(pssm[l])):
                pssm[l][c] = round(math.log2 ((pssm[l][c] + pseudo)/(len(self.seqs) + (len(self.alfabeto)*pseudo)) / (1/len(self.alfabeto))), 2)

        Alfabeto = []
        for z in self.alfabeto:
            Alfabeto.append(z)
        zipped = list(zip(Alfabeto, pssm))
        l = []
        for x in range(1, len(self.seqs[0])+1):
            l.append(x)
        print(f'{"   "}',*l, sep = "    ")
        for x, y in zipped:
            print(x, *y, sep= "  ")

    
    def prob_seq(self, seq):
        """Recebe uma sequência

        Args:
            seq (str): sequência

        Returns:
            int: Probabilidade de gerar uma determinanda sequência
        """
        self.mat = Mat(len(self.alfabeto)+1,len(self.seqs[0])+1)
        prob = 1
        for c in self.seqs:
            for coluna in range(len(self.seqs[0])):
                linha = self.alfabeto.index(c[coluna])
                self.mat[linha + 1][coluna + 1] += 1
                self.mat[0][coluna + 1] = f" {coluna +1} " " " 
                self.mat[0][0] = "." 
                self.mat[linha +1][0] = self.alfabeto[linha]
        for c in range(1,len(self.alfabeto)+ 1):
            for a in range(1,len(self.seqs)+ 1):
                self.mat[c][a] = round(float(self.mat[c][a] + 0.01) / float(len(self.seqs)),3)

        for c in range(1,len(self.seqs[0])+ 1):
            a = self.alfabeto.index(seq[c - 1])
            prob *= self.mat[a+1][c]
        return round(prob,3)


    def seq_mais_provavel(self,pseudo = 0):
        """Recebe pseudo-contagens

        Args:
            pseudo_count (int, optional): _description_. Defaults to 0.

        Returns:
            str: Sequência mais provável
        """
        self.mat = Mat(len(self.alfabeto)+1,len(self.seqs[0])+1)
        resul = ""
        for c in self.seqs:
            for coluna in range(len(self.seqs[0])):
                linha = self.alfabeto.index(c[coluna])
                self.mat[linha + 1][coluna + 1] += 1
                self.mat[0][coluna + 1] = f" {coluna +1} " " " 
                self.mat[0][0] = "." 
                self.mat[linha +1][0] = self.alfabeto[linha] 
        for c in range(1,len(self.alfabeto)+ 1):
            for a in range(1,len(self.seqs)+ 1):
                self.mat[c][a] = round(float(self.mat[c][a] + pseudo + 0.01) / float(len(self.seqs)),3)
        list_temp= []
        list_final = []
        for c in range(1,len(self.seqs)+1):
            for a in range(1,len(self.seqs)+ 1):
                list_temp.append((self.mat[a][c],a))
            list_final.append(max(list_temp))
            list_temp =[]

        for c in list_final:
            for a in range(1,len(self.seqs[0])+1):
                if a == c[1]:
                    resul += self.mat[a][0]
        return resul


import unittest

class Test_Motifs(unittest.TestCase):
    """Classe para a testagem dos métodos 
    """

    def test_seq_mais_provavel(self):
        """Teste para método seq_mais_provavel
        """
        self.assertEqual(Motifs(['ATTG','ATCG','ATTC','ACTC']).seq_mais_provavel(0.5), 'ATTG')


    def test_prob_seq(self):
        """Teste para método prob_seq
        """
        self.assertEqual(Motifs(['ATTG','ATCG','ATTC','ACTC']).prob_seq("ATTG"), 0.284)


unittest.main(argv=[''], exit=False)