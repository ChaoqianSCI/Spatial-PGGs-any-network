# Spatial-PGGs-any-network
This includes a Matlab function for calculating the critical synergy factor of spatial PGGs on any network structure and Python programs for agent-based simulations.

# Files
- everything_r_bc_accu.m - A matlab function. Please input population size N and $N\times N$ adjacency matrix w_ij.
  
  If w_ij(i,j)=1 then i and j are neighbors, and vice versa if w_ij(i,j)=0.

  The function will then return [r_PC,r_DB,r_BD,r_PC_accu,r_DB_accu,r_BD_accu,bc_PC,bc_DB,bc_BD,bc_PC_accu,bc_DB_accu,bc_BD_accu]. "r" means critical synergy factor in spatial PGGs. "bc" means critical benefit-to-cost ratio in spatial DGs. "PC" "DB" "BD" means the three update rules. "accu" means accumulated payoff (if there is no "accu" then it is averaged payoff).

- demo.m - A demonstration for the use of everything_r_bc_accu.m, taking the star graph as an example.

- pgggraph_PC.py - A python program for agent-based simulations of spatial PGGs using averaged payoffs under the PC rule.

- pgggraph_DB.py - A python program for agent-based simulations of spatial PGGs using averaged payoffs under the DB rule.

# System requirements
- We use Matlab 2022a and Python 3.9 in PyCharm.
- We used the mentioned softwares in Windows 11.

# Q&A
For any question about this program, please contact

Chaoqian Wang, Email: CqWang814921147@outlook.com
