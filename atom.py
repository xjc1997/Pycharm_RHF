########updata in 2010-05-11#####################
import periodictable as pt
from definitions import geom
import math
import numpy as np

def do_distmatrix(mol):
    '''本段求取分子中各原子键的距离'''
    mol.dist=np.array([[0.0for i in range (mol.natom)] for j in range (mol.natom)])
    mol.dist2 =np.array([[0.0for i in range (mol.natom)] for j in range (mol.natom)])

    for i in range(mol.natom):
        mol.dist2[i,i]=0.0
        for j in range(i+1,mol.natom):
            mol.dist2[i,j]=sum((mol.x[i,:]-mol.x[j,:])**2)
            mol.dist2[j, i] = mol.dist2[i,j]
            mol.dist[i,j]=math.sqrt(mol.dist2[i,j])
            mol.dist[j,i]=mol.dist[i,j]
    return (mol)
def do_screen(mol):
    nelement_max=54
    tmp=np.array([0 for i in range(nelement_max)])
    ind_atom=np.array([i+1 for i in range(nelement_max)])
    for i in range(mol.natom):
        iatom=int(mol.zatom[i])-1

        tmp[iatom]=tmp[iatom]+1
    mark=(tmp>0)
    mark=mark.tolist()

    mol.ntype=mark.count(True)                    #统计原子数目
    mol.ztype=[0.0 for i in range(mol.ntype)]
    mol.itype_number = [0for i in range(mol.ntype)]
    mol.iatom_to_type = [0 for i in range(mol.natom)]

    mol.ztype=ind_atom[mark]              #由bool值输出原子种类（元素表示核数）
    mol.itype_number=tmp[mark]            #有bool值给出原子个数
    for i in range(mol.ntype):
        intzatom=list(map(int,mol.zatom))
        for j in range(mol.natom):
            if(intzatom[j]==(mol.ztype[i])):
                mol.iatom_to_type[j] = i                     #s输入分子中各个原子在type里的序列
    return(mol)
def print_atoms(mol):                                               #用于打印分子中的原子
    print("NAtoms :",mol.natom)
    print("Natomic type :",mol.ntype)
    print("Element","Number")
    for i in range(mol.ntype):
        name=pt.elements[mol.ztype[i]].symbol
        print(name,'     ',mol.itype_number[i])

def readatoms(filename,mol,print_):
    try:
        f = open(filename,'r+')
        f.close()
    except IOError:
        print("input File is not accessible.")
    ###以上为检查file文件是否存在#############################
    bohr_A=0.52917721092
    mol=geom()
    f = open(filename, 'r+')
    N=int(f.readline())
    f.readline()
    mol.natom=N
    mol.zatom=np.array([0.0 for i in range(N) ])
    mol.x=np.array([[0.0 for i in range(3)] for i in range(N)])  #0.0表示浮点数
    for i in range(N):
        a=f.readline().strip('\n').split()                         #去除每一行的换行符与空格
        symbol=a[0]                                                #读取原子符号名
        a=np.array(a[1:])                                          #去除原子符号名
        a=[float(num) for num in a ]
        mol.x[i,:]=a
        symbol=pt.elements.symbol(symbol)
        mol.zatom[i]=symbol.number
    mol.x=mol.x/bohr_A
    mol=do_distmatrix(mol)
    #print(mol.dist2)                                        #测试各原子间输出距离4
    mol=do_screen(mol)
    print(mol.iatom_to_type)
    f.close()
    if(print_==1):                                         #打印分子中的原子
        print("==============atom list ===============")
        print_atoms(mol)
    return(mol)




