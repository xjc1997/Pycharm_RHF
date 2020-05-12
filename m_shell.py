from atom import geom
import sys
import periodictable as pt
import math
import numpy as np
import os
class shell:
    def __init__(self):
        self.am=0
        self.bool_sp=True
        self.itype=0
        self.ncontr=0
        self.alpha=[]
        self.coeff=[[]]
class shell_pair:
    def __init__(self):
        self.o2a=[]
        self.o2b=[]
        self.zeta=[]
        self.oozeta=[]
        self.oo2zeta=[]
        self.rho=[]
        self.boab=[]
        self.aoab=[]

class shell_set:
    def __init__(self):
        self.basis_set_name=""
        self.nshell=0
        self.nbf_shell=[]
        self.nbf_atype=[]
        self.ishell_itype=[]
        self.nshell_atype=[]
        self.shells=[shell()]
        self.shell_pair_data=[[shell_pair()]]


def number_basis_function_am(guassiantype,am):
    if(guassiantype=='CART'):
        #print(guassiantype)                                    #输入基组类型
        number_basis_function_amz=(am+1)*(am+2)/2
    return(number_basis_function_amz)


def doublefactorial(n):
    if n <= 0:
        return 1
    else:
        return n * doublefactorial(n - 2)

def renorm(coeff,am,alpha):
    sqrt_Pi_cubed = 5.568327996831707
    # if(any(np.array(alpha) < 0.0)):
    #      print('coefficent is error ,please check it')
    #      sys.exit()
    alpha=np.array(alpha)

    tmp=2.0*alpha
    tmp=tmp**(am+1.5)
    double_factorial=doublefactorial(2*am-1)
    #print(double_factorial)
    tmp=np.sqrt((2.0**am*tmp)/((sqrt_Pi_cubed)*double_factorial))
    norm_coeff=tmp*coeff
    return (norm_coeff)
#############################################################################
def init_shell(ncontr,am,itype,alpha,coeff1,coeff2):
    "本段的作用是将输入文件的内容转入到SHELL类中的实例中去"
    sh=shell()            #创建一个实例
    if(am<0):
        print('DATA ERROR IN SHELL')
        sys.exit()
    sh.alpha=[0.0 for i in range(ncontr)]
    if(am==10):
        sh.coeff=np.array([[0.0 for i in range(ncontr) ]for j in range(2)])
        sh.bool_sp=True
    else:
        sh.coeff=np.array([0.0for i in range(ncontr)])
        sh.bool_sp=False
    sh.ncontr=ncontr
    sh.am=am
    sh.itype=itype
    sh.alpha[:]=alpha[:]
    if(am!=10):
        sh.coeff=renorm(coeff1,am,alpha)
    else:
        sh.coeff[0,:]=renorm(coeff1,0,alpha)
        sh.coeff[1, :] = renorm(coeff1, 1, alpha)
        sh.am=1
    return(sh)
########################################## #########################




def read_shell_set(basis_path,basis_name,mol,basis,print_):
    basis.nshell=0
    for itype in range(mol.ntype):
        symbol=pt.elements[mol.ztype[itype]].symbol
        basis_filename=basis_path+"/"+symbol+"_"+basis_name                                               ##确定基组所在路径
        try:
            f = open(basis_filename, 'r+')
            f.close()
        except IOError:
            print("basis input File"+basis_filename+" is not accessible.")     ##基组文件debug

        #############以上为检查文件是否为空#######################################
        f = open(basis_filename, 'r+')
        nbf_file=int(f.readline())                ##读取第一行以获取壳层数
        if(nbf_file<1):
            print(basis_filename+' is free,please check it')
            sys.exit()
        basis.nshell=basis.nshell+nbf_file
        f.close()
      ####以上为获取应设置数组的大小#################
    basis.shells=[shell() for i in range(basis.nshell)]
    basis.nbf_shell=[0 for i in range(basis.nshell)]
    basis.nbf_atype=[0 for i in range(mol.ntype)]
    basis.ishell_atype = [0 for i in range(mol.ntype)]
    basis.nshell_atype = [0 for i in range(mol.ntype)]

    ishell=1
    for itype in range(mol.ntype):
        symbol = pt.elements[mol.ztype[itype]].symbol
        basis_filename = basis_path + "/" + symbol + "_" + basis_name
        f = open(basis_filename, 'r+')
        nbf_file=int(f.readline())
        ishell_l=ishell
        basis.ishell_atype[itype]=ishell_l
        for ibf_file in range(nbf_file):
            ng,am_tmp=f.readline().split()
            ng=int(ng)
            am_tmp=int(am_tmp)
            if (ng<1 ):
                print('error happen in basis set file'+basis_filename+'please check it')
                sys.exit()
            ####################以上为检查输入文件#########################################
            alpha=[0.0 for i in range(ng)]
            coeff1 = [0.0 for i in range(ng)]
            coeff2 = [0.0 for i in range(ng)]
            if(am_tmp<10):
                for ig in range(ng-1,-1,-1):
                    alpha[ig],coeff1[ig]=f.readline().split()
                    alpha[ig]=float(alpha[ig])
                    coeff1[ig]=float(coeff1[ig])
            else:
                for ig in range(ng-1,-1,-1):
                    alpha[ig],coeff1[ig],coeff2[ig]=f.readline().split()
                    alpha[ig] = float(alpha[ig])
                    coeff1[ig] = float(coeff1[ig])
                    coeff2[ig] = float(coeff2[ig])
            #print('basis shell number',ishell)
            basis.nbf_shell[ishell-1]=number_basis_function_am('CART',am_tmp)
            basis.shells[ishell-1]=init_shell(ng,am_tmp,itype,alpha,coeff1,coeff2)
            ishell=ishell+1
        basis.nshell_atype[itype]=ishell-ishell_l
        basis.nbf_atype[itype]=sum(basis.nbf_shell[ishell_l-1:ishell-1])
        #print(itype,ishell_l-1,ishell-1,basis.nbf_shell,basis.nbf_atype[itype])
        f.close()
        return (basis)


