'''
script for creating commands
python3 survival_auto.py abc xyz
where
abc: x-position (0.8 for now)
xyz: alpha angle (from 0 to 180)
output should be to auto.sh 
'''
for i in range(1,10):
    for j in range(1,180):
        print('python3 asymmetric.py '+str(i/10)+' '+str(j)+' >> resultados.txt' )

