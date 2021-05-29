# Contributing to SpectroscoPyx

There are numerous ways to contribute to SpectroscoPyx, including by
providing code and documentation, suggesting and discussing ideas,
submitting issues and bug reports, and engaging the broader plasma
physics community.  

**Imposter syndrome disclaimer**: We want your help. No, really.

There may be a little voice inside your head that is telling you that you're not
ready to be an open source contributor; that your skills aren't nearly good
enough to contribute. What could you possibly offer a project like this one?

We assure you - the little voice in your head is wrong. If you can write code at
all, you can contribute code to open source. Contributing to open source
projects is a fantastic way to advance one's coding skills. Writing perfect code
isn't the measure of a good developer (that would disqualify all of us!); it's
trying to create something, making mistakes, and learning from those
mistakes. That's how we all improve, and we are happy to help others learn.

Being an open source contributor doesn't just mean writing code, either. You can
help out by writing documentation, tests, or even giving feedback about the
project (and yes - that includes giving feedback about the contribution
process). Some of these contributions may be the most valuable to the project as
a whole, because you're coming to the project with fresh eyes, so you can see
the errors and assumptions that seasoned contributors have glossed over.

*This disclaimer was originally written by
[Adrienne Lowe](https://github.com/adriennefriend) for a
[PyCon talk](https://www.youtube.com/watch?v=6Uj746j9Heo), and was adapted by 
[yt](https://github.com/yt-project/yt) in their README file based on its use 
in the README file for the [MetPy project](https://github.com/Unidata/MetPy).
It was then adapted by SpectroscoPyx.*

## Sharing ideas

There are several methods of communication that are being used in the
early stages of SpectroscoPyx development:

* [Signing up for the PlasmaPy email
  list](https://groups.google.com/forum/#!forum/plasmapy) will allow
  you to participate in broader discussions and keep up with the
  latest developments.

* The [SpectroscoPyx repository on
  GitHub](https://github.com/PlasmaPy/SpectroscoPyx) is the best place to
  [submit issues](https://github.com/PlasmaPy/SpectroscoPyx/issues) and
  review [pull requests](https://github.com/PlasmaPy/SpectroscoPyx/pulls).

* The PlasmaPy [Matrix](https://riot.im/app/#/room/#plasmapy:matrix.org) or 
  [Gitter](https://gitter.im/PlasmaPy/Lobby) joint channel
  is a great place to have informal conversations, coordinate efforts,
  and share ideas.  
* We have biweekly telecons which are announced on the email list.

## Contributing code or documentation to SpectroscoPyx

If you see something you'd like to work on amongst our
[issues](https://github.com/PlasmaPy/SpectroscoPyx/issues), start hacking away on
 that! However, please announce your intent first in the relevant issue to 
 make sure there is no work duplication.
 
Please note that SpectroscoPyx has a [Code of Conduct](CODE_OF_CONDUCT.md).

Issues marked by the community as *help wanted* mean just that - either they're good contributions for outsiders or there's an issue in the ongoing work that requires a second opinion. Please consider these first!

### Preliminaries

Work on SpectroscoPyx is done via GitHub, so you'll need a
[(free) account](https://github.com/join?source=header-home).
If you are new to [git](https://git-scm.com/), helpful resources include
documentation on [git
basics](https://git-scm.com/book/en/v2/Getting-Started-Git-Basics) and
an [interactive git
tutorial](https://try.github.io/levels/1/challenges/1).  You must also
[install
git](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git)
locally on your computer.  We highly recommend getting familiar with
git by going through these tutorials or a [Software
Carpentry](https://software-carpentry.org/) workshop prior to making
code contributions.

### Virtual Environments

Before you grab SpectroscoPyx from GitHub, you are going to want to setup a sensible
directory structure and a virtual environment. The virtual environment will
allow you to import, run, and test your development version of SpectroscoPyx without
contaminating or conflicting with other version of SpectroscoPyx or other packages
that may be on your system.

If you are running the Anaconda suite and want to use virtualenv to setup your
virtual environment, you will have to let the system know where the Python
interpreter can be found. On Linux this is done with (for example, assuming
having installed Anaconda into `~/anaconda3`):

```ShellSession
export LD_LIBRARY_PATH="$HOME/anaconda3/lib/"
```

Exporting the library path to the dynamic linker will only last for the
duration of the current shell session.

Next you should create a sensible directory structure. Something like:

```ShellSession
mkdir ~/programming/spectroscopyx/
```

You need to make sure that the directory path names don't contain any spaces,
otherwise virtualenv will throw an error. Now move into the directory and
create the virtual environment

```ShellSession
cd ~/programming/spectroscopyx/
virtualenv -p python3 .
```

Your virtual environment should now be created. If you run `ls` you will notice
that virtualenv has created a number of subdirectories: `bin/`, `lib/`, and
`include/`. To activate the virtualenv you will run:

```ShellSession
source ./bin/activate
```

You will have to add the python library directory to LD_LIBRARY_PATH, as
described in a previous step, prior to activating the virtualenv for every new
shell session.  You should now see that your shell session is prepended with
(spectroscopyx), like so:

```ShellSession
(spectroscopyx) user@name:~/programming/spectroscopyx$ 
```

This indicates that the virtualenv is running. Congratulations!
When your're done working on SpectroscoPyx, you can deactivate the virtualenv by
running

```ShellSession
source deactivate
```

If you are running virtualenv, then in the next step you will want to clone
SpectroscoPyx while in the `~/programming/spectroscopyx` directory. This will create a
subdirectory `~/programming/spectroscopyx/SpectroscoPyx/` which will prevent the package
from being contaminated will all those `lib/`, `bin/`, `include/` directories
which virtualenv generated. Alternatively, you can setup the clone first, and
then setup the virtualenv inside `SpectroscoPyx/`, making sure to add those
virtualenv directories into a .gitignore file so that they don't get pushed
upstream.

### Forking and cloning SpectroscoPyx

After creating your GitHub account, go to the [main
repository](https://github.com/PlasmaPy/SpectroscoPyx) and **fork a copy of
SpectroscoPyx to your account**.

Next you must **clone your fork to your computer**.  Go to the
directory that will host your SpectroscoPyx directory, and run one of the
following commands (after changing *your-username* to your username).
If you would like to use HTTPS (which is the default and easier to set
up), then run:

```ShellSession
git clone https://github.com/your-username/SpectroscoPyx.git
```

SSH is a more secure option, but requires you to [set up an SSH
key](https://help.github.com/articles/generating-a-new-ssh-key-and-adding-it-to-the-ssh-agent/)
beforehand.  The equivalent SSH command is:

```ShellSession
git clone git@github.com:your-username/SpectroscoPyx.git
```

After cloning, we must tell git where the development version of
SpectroscoPyx is by running:

```ShellSession
git remote add upstream git://github.com/PlasmaPy/SpectroscoPyx.git
```

To check on which remotes exist, run `git remote -v`.  You should get
something like this:

```ShellSession
origin		git@github.com:namurphy/SpectroscoPyx.git (fetch)
origin		git@github.com:namurphy/SpectroscoPyx.git (push)
upstream	git@github.com:PlasmaPy/SpectroscoPyx.git (fetch)
upstream	git@github.com:PlasmaPy/SpectroscoPyx.git (push)
```

### Setting up SpectroscoPyx for testing

Now that you have SpectroscoPyx on your local computer and you have a virtual
environment, you will want to "install" this development version of SpectroscoPyx
along with its dependencies. Start by activating your virtual environment. Next
you want install the SpectroscoPyx dependencies. One way to do this is to do

```ShellSession
(spectroscopyx) user@name:~/programming/spectroscopyx$ pip install -r requirements/base.txt
```

Next, setup the development version of SpectroscoPyx which you just cloned by moving
into the root directory of the cloned repo and running the setup.py script
there:

```ShellSession
(spectroscopyx) user@name:~/programming/spectroscopyx/SpectroscoPyx$ python setup.py develop
```

You should now be all set to run development versions of SpectroscoPyx modules via
`import SpectroscoPyx` in your test scripts!


### Branches, commits, and pull requests

Before making any changes, it is prudent to update your local
repository with the most recent changes from the development
repository:

```ShellSession
git fetch upstream
```

Changes to SpectroscoPyx should be made using branches.  It is usually best
to avoid making changes on your main branch so that it can be kept
consistent with the upstream repository.  Instead we can create a new
branch for the specific feature that you would like to work on:

```ShellSession
git branch *your-new-feature*
``` 

Descriptive branch names such as `grad-shafranov` or
`adding-eigenfunction-poetry` are helpful, while vague names like
`edits` are considered harmful.  After creating your branch locally,
let your fork of SpectroscoPyx know about it by running:

```ShellSession
git push --set-upstream origin *your-new-feature*
``` 

It is also useful to configure git so that only the branch you are
working on gets pushed to GitHub:

```ShellSession
git config --global push.default simple
```

Once you have set up your fork and created a branch, you are ready to
make edits to SpectroscoPyx.  Switch to your new branch by running:

```ShellSession
git checkout *your-new-feature*
```

Go ahead and modify files with your favorite text editor.  Be sure to
include tests and documentation with any new functionality.  We also
recommend reading about [best practices for scientific
computing](https://doi.org/10.1371/journal.pbio.1001745).  SpectroscoPyx
uses the [PEP 8 style guide for Python
code](https://www.python.org/dev/peps/pep-0008/) and the [numpydoc
format for
docstrings](https://github.com/numpy/numpy/blob/main/doc/HOWTO_DOCUMENT.rst.txt)
to maintain consistency and readability.  New contributors should not 
worry too much about precisely matching these styles when first 
submitting a pull request, as the [PEP8 Speaks](http://pep8speaks.com/)
GitHub integration will check pull requests for PEP 8 compatibility, and
further changes to the style can be suggested during code review.

You may periodically commit changes to your branch by running

```ShellSession
git add filename.py
git commit -m "*brief description of changes*"
```

Committed changes may be pushed to the corresponding branch on your
GitHub fork of SpectroscoPyx using 

```ShellSession
git push origin *your-new-feature* 
```

or, more simply,

```ShellSession
git push
```

Once you have completed your changes and pushed them to the branch on
GitHub, you are ready to make a pull request.  Go to your fork of
SpectroscoPyx in GitHub.  Select "Compare and pull request".  Add a
descriptive title and some details about your changes.  Then select
"Create pull request".  Other contributors will then have a chance to
review the code and offer contructive suggestions.  You can continue
to edit the pull request by changing the corresponding branch on your
SpectroscoPyx fork on GitHub.  After a pull request is merged into the
code, you may delete the branch you created for that pull request.
