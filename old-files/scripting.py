from subprocess import call, Popen, PIPE


call(['sbatch', 'shell-executer.sh'])

#cmd = "FRED=HELLO"

#p = Popen(cmd, stdin=PIPE, stdout=PIPE)
#stdout = p.communicate()[0]
#print stdout
