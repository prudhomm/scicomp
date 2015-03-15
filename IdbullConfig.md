# Introduction #

idbull is a ia64 arch with 8 processors. The system is a modified redhat linux distribution.



# Login info #

Perform the following operations:

  * download the .toprc from the main project page to your $HOME directory for a proper 'top' utility configuration

  * create the file `$HOME/.bash_profile`  if not already created and make sure it contains the following
```
# .bash_profile

# Get the aliases and functions
if [ -f ~/.bashrc ]; then
        . ~/.bashrc
fi

# User specific environment and startup programs

PATH=/usr/ljk/bin:/usr/local/mpich/bin:$PATH:$HOME/bin

export PATH
unset USERNAME

LD_LIBRARY_PATH=/usr/ljk/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH
```

  * at the moment you can only connect from the UFR IMA (UFR Phys access should be opened soon).  First log on `mandelbrot.e.ujf-grenoble.fr` then on `idbull.ujf-grenoble.fr`

## ssh connection ##

to ease the connection on idbull it is advised to use the ssh keys

  * create the ssh keys (and enter your passphrase)
```
 ssh-keygen -t dsa
```

  * copy .ssh/id\_dsa.pub on idbull:~/.ssh/authorized\_keys

  * type
```
ssh-add
```
and type your passphrase

  * you do that on the machines you are logged on.

  * then no need to retype the passphrase anymore