---
title: Installing Magic
author: Tim Stuart
date: '2017-03-03'
comments: true
---

# Installing magic

Recently a method for imputing single cell gene expression matricies
was posted
on [biorxiv](http://biorxiv.org/content/early/2017/02/25/111591) by
David van Dijk et al., called magic (Markov Affinity-based Graph
Imputation of Cells). I've been analysing single cell RNA-seq data
recently, and this method looks like it could be useful when trying to
find co-transcriptional networks, as single cell data suffers from
dropout which makes finding co-transcriptional networks hard.

I had lots of problems getting magic installed and running, so will
document them here for future reference.

Firstly, there seems to be an error in the `setup.py` script where it
looks for a non-existent `/data` directory that should contain test
data. Running `pip3 install .` as instructed in the readme resulted in
an error, and I have raised
an [issue](https://github.com/pkathail/magic/issues/12) on
github. Commenting out the last few lines of the `setup.py` script
seemed to provide a temporary fix.

The next problem was getting `Tk` to work properly with python. `Tk`
is a GUI toolkit and not park of python itself. Chances are
that `Tk` is installed somewhere on your computer, and the problem is that
python doesn't know where it is. After trying lots of different
things, the solution I found was to install python3 using the mac
installer and launching IDLE, as this finds and links the `Tk`
installation with python at runtime. From the
python [website](https://www.python.org/download/mac/tcltk/):

> The Python for Mac OS X installers downloaded from this website
> dynamically link at runtime to Tcl/Tk macOS frameworks. 

I then found the path of the newly installed python3 (it was symlinked
to `/usr/local/bin/python3` for me) and used this to create a new
virtualenv:

```
mkvirtualenv -p /usr/local/bin/python3 py3
```

I then installed magic again in the virtualenv (from the github repo):

```
pip3 install .
```

Next I installed all the jupyter stuff. It's important to link the
right ipython kernel to the jupyter notebook, otherwise it will seem
like you still don't have access to `Tk`, even though at this point
you can sucessfully `import tkinter` in python3. To do this, install
jupyter, ipython, the ipython kernel, and then link the
kernel:

```
pip install jupyter ipython ipykernel
python3 -m ipykernel install --user
```

Now you should be able to import magic without any problems, and use
it in a jupyter notebook. You can also start the GUI by running the `magic_gui.py` script:

```
python3 magic_gui.py
```

<img src="/assets/magic.png" alt="Drawing" style="width: 80%;"/>
