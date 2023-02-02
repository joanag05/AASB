import io

class Mat:
  """Implementação de uma matriz
  """
  def __init__(self, rows, cols):                                         
      """Construtor da classe

      Args:
          rows (int): _linhas da matrz
          cols (int): _colunas da matriz
      """
      self.mat = [[0 for c in range(cols)]
                    for r in range(rows)]

  def numRows (self): return len(self.mat)

  def numCols (self): return len(self.mat[0])

  def __str__(self):
        """"Devolve a matriz como uma string
        """
        return '\n'.join(' '.join(str(val) for val in row)
                         for row in self.mat)
    
  def __repr__(self):
        return str(self)
    
  def __getitem__ (self, n):
        """Interface para a indexação []
        """
        return self.mat[n]


def subst(x,y):
  """Função que permite substituir valores de match e mismatch no alinhamento local

  Args:
    x (str): elemento da s1
    y (str): elemento da s2 
  Returns:
    int: valores de match(+4) e mismatch(-2)
    """  
                                                        
  return 4 if x == y else -2                                          

class NW:
  """Alinhamento global entre duas sequências usando valores de match e mismatch pré defenidos

  """
  def __init__(self, s1, s2, g = -2):
    """Contrutor da classe NW

    Args:
        s1 (str): 1ª sequência
        s2 (str): 2ª sequência
        g (int): Gap Penalty. Default = -2
    """
    self.s1 = s1
    self.s2 = s2
    self.mat = Mat(len(s1) + 1, len(s2) + 1)
    self.tr  = Mat(len(s1) + 1, len(s2) + 1)
    
    for L in range(len(s1)):
      self.mat[L + 1][0] = g * (L + 1)
      self.tr[L + 1][0]  = 'C'

    for C in range(len(s2)):
      self.mat[0][C + 1] = g * (C + 1)
      self.tr[0][C + 1]  = 'E'

    for L, x1 in enumerate(s1):
      for C, x2 in enumerate(s2):
        possiveis = [
            self.mat[L  ][C    ] + subst(x1, x2),   
            self.mat[L+1][C    ] + g,              
            self.mat[L  ][C + 1] + g,               
        ]
        dirs = "DEC"

        self.mat[L + 1][C + 1] = max(possiveis)
        self.tr[L + 1][C + 1] = dirs[possiveis.index(self.mat[L + 1][C + 1])]
  
  def align_score(self):
    """ Função que recebe os alinhamentos

    Retorna:
        int: O score máximo do alinhamento NW entre duas sequências 
    """
    return self.mat[len(self.s1)][len(self.s2)]

  def rebuild(self):
    """Função que implementa o processo de reconstrução e recuperação do 
    alinhamento ótimo entre duas sequências

    Retorna:
        str: Alinhamento Ótimo entre as duas sequências
    """
                                              
    L = len(self.s1)
    C = len(self.s2)
    S1 = ""
    S2 = ""
    
    dirs = {
        'D' : (-1, -1),               
        'E' : ( 0, -1),               
        'C' : (-1,  0)               
    }

    while L > 0 or C > 0:
      DL, DC = dirs[self.tr[L][C]]


      if self.tr[L][C] == "D":
        S1 = self.s1[L - 1] + S1
        S2 = self.s2[C - 1] + S2
      elif self.tr[L][C] == "E":
        S1 = '-' + S1
        S2 = self.s2[C - 1] + S2
      else:
        S1 = self.s1[L - 1] + S1
        S2 = '-' + S2        

      L += DL
      C += DC

    return S1, S2

  def __repr__(self):
    """Função que premite ver a representação das matrizes S e T

    Retorna:
        _type_: Matriz S e Matriz T
    """
    cols = "-" + self.s2
    lins = "-" + self.s1
    with io.StringIO("") as S:
      print("\n Matriz S :")
      print(' ', *cols, sep = '   ', file = S)
      for L, linha in zip(lins, self.mat):
        print(L, *[f'{x:3d}' for x in linha], file = S)

      print('\n Matriz T :', file = S)

      print(' ', *cols, file = S)
      for L, linha in zip(lins, self.tr):
        print(L, *linha, file = S)

      return S.getvalue()



import unittest

class Test_NW_ALIGNMENT(unittest.TestCase):
    """Classe que permite testagem do código do alinhamento NW
    """
    def test_score(self):
        """Testar a função score
        """
        self.assertEqual(NW("ACHDHGT","THHDACCGT").align_score(), 6)
        self
    
    def test_align_rebuild(self):
        """Testar a função rebuild
        """
        self.assertEqual (NW("ACHDHGT","THHDACCGT").rebuild(), ('ACHD--HGT', 'THHDACCGT'))




    
if __name__ == '__main__':
    unittest.main()