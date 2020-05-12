#updata : 2010-05-10  xjc1997@mail.ustc.edu.cn
import math
Ha_eV=27.21138505
bohr_A=0.52917721092
c_speedlight=137.035999074
pi=3.14159265358979323
pi2=pi**2
sqrtpi=math.sqrt(pi)
sqrt2=math.sqrt(2)
im=0.0+1.0j
npar_arth = int(16)
npara2_arth = int(8)
npar_poly = int(8)


class geom:
    def __init__(self):
        self.natom=0
        self.ntype=0
        self.zatom=[]
        self.ztype=[]
        self.x=[[]]
        self.dist=[[]]
        self.dist2=[[]]
        self.iatom_to_type=[]
        self.itype_number=[]

