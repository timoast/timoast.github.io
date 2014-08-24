# Useful Python Code

### Take command-line arguments
```python
from sys import argv


def checkArgs(arg1, arg2):
    """
    arg1 is short arg, eg h
    arg2 is long arg, eg host
    """
    args = argv[1:]
    if arg1 in args:
        index = args.index(arg1)+1
        variable = args[index]
        return variable
    elif arg2 in args:
        index = args.index(arg2)+1
        variable = args[index]
        return variable
    else:
        variable = raw_input("\nEnter {arg2}: ".format(arg2=arg2))
        return variable

```

## List Comprehension
If ever declaring and empty list then using `.append`, can be done with list comprehension.
```python
list3 = []
for x in list2:
    if x[0] in list1:
        list3.append(x)

list3 = [x for x in list2 if x[0] in list1]
```

## Printing
### Print in colour
```python
class colour:
    PURPLE = '\033[95m'
    CYAN = '\033[96m'
    DARKCYAN = '\033[36m'
    BLUE = '\033[94m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    END = '\033[0m'

print colour.CYAN + 'This is blue' + colour.END
```

### Format text wrapping correctly in terminal
```python
import textwrap
print textwrap.dedent('print message')
```

## Plotly
### Import statements
```python
import plotly.plotly as py
from plotly.graph_objs import *
import plotly.tools as tls
```

### Create credentials file
```python
import plotly.tools as tls
tls.set_credentials_file(username="your_username", api_key="your_api_key")
```

### Sign in with credentials file
```python
my_creds = tls.get_credentials_file()
py.sign_in(my_creds['username'], my_creds['api_key'])
```

### Layout settings
```python
l = Layout(title= 'Title',
  yaxis= YAxis(title= 'X-axis title'),
  xaxis= XAxis(title= 'Y-axis title'),
  barmode= 'overlay', # group, stack
  showlegend= False,
  width= 500,
  height= 500,
  autosize= False,
  margin= Margin(
  l= 50,
  r= 50,
  b= 100,
  t= 100,
  pad= 4
  ))
```

### Figures
```python
dataset1 = Histogram(x= dataset1, name= 'Data', opacity= 0.75)
dataset1 = Histogram(x= dataset2, name= 'Data', opacity= 0.75)
data = Data([dataset1, dataset2])
fig = Figure(data= data, layout= l)
py.iplot(fig) # interactive, inline
plot_url = py.plot(fig, auto_open= False) # returns plot url
```

## File reading / writing
### Read lines in tsv/csv file
```python
with open('merged_dmrs.bed', 'r') as f:
    for line in f:
        line = line.rsplit() # each item in line is now part of a list
```

## Mapping
Functions for mapping gDNA reads
```python
import os
from subprocess import call


def peMap(proc):
    """
    paired-end reads
    """
    for filename in os.listdir('.'):
        if filename.endswith("_1.fastq"):
            sam = filename.split('_')
            sam = sam[0]
            call(["bowtie2", "-p{x}".format(x=proc), "-q", "-x /dd_stage/userdata/lister/data/genomes/bowtie2_indexes/tair9",
                  "-1 {f}".format(f=sam + '_1.fastq'), "-2 {f}".format(f=sam + '_2.fastq'), "-S {s}".format(s=sam + '.sam')])
        else:
            pass


def seMap(proc):
    """
    single-end reads
    """
    for filename in os.listdir('.'):
        if filename.endswith(".fastq"):
            sam = filename.split('.')
            sam = sam[0]
            call(["bowtie2", "-p{x}".format(x=proc), "-q", "-x /dd_stage/userdata/lister/data/genomes/bowtie2_indexes/tair9",
                  "-U {f}".format(f=filename), "-S {s}".format(s=sam + '.sam')])
        else:
            pass
```

## Split sra files into fastq
```python
import os
from subprocess import call


def fastqSplit():
    for filename in os.listdir('.'):
        if filename.endswith('.sra'):
            print 'processing {n}'.format(n=filename)
            call(['fastq-dump', '--split-3', '-v', filename])
        else:
            pass
```
