from tabulate import tabulate
from texttable import Texttable
import latextable

def build_La_Tex_table(rows, caption = "A comparison of experiment results with different values of $k$"):
    
    table = Texttable()
    table.set_cols_align(["c"] * 4)
    table.set_deco(Texttable.HEADER | Texttable.VLINES)
    table.add_rows(rows)


    print('\nTexttable Latex:')
    print(latextable.draw_latex(table, caption=caption))