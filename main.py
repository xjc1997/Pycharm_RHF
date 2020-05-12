from atom import readatoms
from definitions import geom
from m_shell import read_shell_set
from m_shell import shell_set
from m_shell import init_basis_set
print_=1                      #打印参数
filename="out.txt"
outfile=open(filename,'w+')
geom_path="ch4.txt"
mol=geom()

mol=readatoms(geom_path,mol,print_)
basis_path="."
basis_name="STO-3G"
shells=shell_set()
shells=read_shell_set(basis_path,basis_name,mol,shells,print_)
init_basis_set(mol,shells,basis,print_)