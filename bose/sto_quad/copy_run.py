import os
for i in range(10):
    cmd = r'cd ' + str(i) + '&&'
    cmd += r'echo new > out' + '&&'
    cmd += r'cp ../decompose_any_spe.py ./&&'
    cmd += r'cp ../r.sh ./&&'
    cmd += r'python3 decompose_any_spe.py >> out&'

    os.system(cmd)
