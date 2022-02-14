# Guidance on how to contribute

> All contributions to this project will be released to the public domain.
> By submitting a pull request or filing a bug, issue, or
> feature request, you are agreeing to comply with this waiver of copyright interest.
> Details can be found in our [TERMS](TERMS.md) and [LICENSE](LICENSE).

There are two primary ways to help:
 - Using the issue tracker, and
 - Changing the code-base.

## Using the issue tracker
<!--- TODO: Change this to refer specifically to the GitHub issues feature, if permitted under organization policy-->
Use the issue tracker to suggest feature requests, report bugs, and ask questions.
This is also a great way to connect with the developers of the project as well
as others who are interested in this solution.

Use the issue tracker to find ways to contribute. Find a bug or a feature, mention in
the issue that you will take on that effort, then follow the _Changing the code-base_
guidance below.

## Changing the code-base

Generally speaking, you should fork this repository, make a 
new branch, make changes in your own fork/branch, and then 
submit a pull request to incorporate changes into the main 
codebase. All new code *should* have associated 
unit tests that validate implemented features and the presence 
or lack of defects. In almost all cases, the submitted code 
should have some kind of demonstration notebook
to go with it so that we can help others experiment with our project.

Additionally, the code should follow any stylistic and 
architectural guidelines prescribed by the project. 
We request that you use ['Black'](https://pypi.org/project/black/)
before submitting any pull requests. 


<!--- TODO: Consider using or merging this GitHub description with https://github.com/NOAA-OWP/DMOD/blob/master/doc/GIT_USAGE.md --->
## The process in summary, i.e., TL;DR
* Fork the t-route repo on GitHub and clone the new fork  to your development environment

```
git clone github.com/<githubusername>/t-route.git
```

* Make a new branch and make your changes to the source files in the branch

* Add changed files to the commit ... `git add <list of your changed files>`

* ... and commit them: `git commit`

* Push the accumulated commits to your fork/branch in GitHub: `git push`

* Open GitHub and issue a pull request.

* **IMPORTANT** Keep your master branch code up-to-date with the upstream repository using the following commands. As commits from other developers accumulate in the master, merge these into your fork prior to issuing a pull request.

```
git remote add upstream https://github.com/NOAA-OWP/t-route.git
git fetch upstream && git pull --rebase upstream master

```
## A more detailed step-by-step guide to contributing via GitHub

### 1. On GitHub, create a fork and clone it to your local directory 
A step-by-step guide illustrating how to fork and clone 
a repository can be found 
[here](https://help.github.com/en/github/getting-started-with-github/fork-a-repo). 
Once a new fork is created, clone it to an empty local directory 
of your choice (if the directory does not exist already, 
`git clone` will create it):

```
# navigate to a local directory where you want to store the repo
cd <local/directory/of/your/choice>

# clone your forked repo
git clone github.com/<githubusername>/t-route.git
 ```  
____

### 2. Set up the Git configuration parameters

Set up your user name and email in your local repo config

```
git config user.name "Your GitHub username"
git config user.email "Your email address"
```
Alternatively, you can change your machine global git configs: 
    
```
git config --global user.name "Your Github username"
git config --global user.email "Your email address"
```
____

<!--- Please check if this is the default behavior. If so, 
we can probably remove this step --->
### 3. Make your fork a remote 
It is common practice to name the remote, `origin`

```
git remote add origin master https://github.com/<user.name>/t-route.git
```
____

### 4. Add the upstream t-route repo as a second remote

To get the latest updates from the main t-route codebase, 
we need to add it as a second remote. It is common to 
name this remote, `upstream` 

```
git remote add upstream https://github.com/NOAA-OWP/t-route.git
```

____


### 5. Incorporate upstream changes into your fork

To incorporate any changes from the main codebase into your
local repo, `checkout` your `master` branch and do a 
fetch-and-rebase using the `pull` command.

```
# checkout the master branch
git checkout master

# fetch and rebase changes 
git pull --rebase upstream master
```

Finally, use the `push` command to build these changes into 
your fork.

```
# assuming you still have `master` checked out
git push
```

The above fetch-and-rebase process is recommended prior to 
creating a new development branch. This will ensure that your 
new additions are relevant to the latest version of the codebase. 

____


### 6. Create a new development branch

Branches are parallel versions of the repository that allow you 
to make changes and tinker as much as you would like without 
disrupting the master codebase. It is recommended to create a 
new branch for each bite-sized contribution you make or issue 
you respond to. Here is how to create a new branch and 
automatically switch to it:

```
git checkout -b <new-branch-name>
```

The above command will create *and* checkout the new branch. 
Otherwise, the `checkout` command is used to navigate between 
branches. It updates files in your working directory to match 
the version stored in that branch. 

```
git checkout <branch-name>
```
____

### 7. Incorporate upstream changes into your development branch

Any time you update `master` with changes from the `upstream` 
remote (step 5), it is advised to rebase any local development
branches, too. Rebasing ensures that your changes are based 
from what everyone else has already done. 

```
# check out the branch you would like to rebase
git checkout <branch-name>

# rebase the branch 
git rebase master
```
____

### 8. Tinker in your development branches. 
In your local development branches, explore the code, make changes, 
and tinker at will. **Be creative!**

____

### 9. Add new and changed files to the staging area

```
git add <list of your changed files>
```
____

### 10. Commit your changes to the checked out development branch
Commits should be accompanied by a concise and brief comment. 
Like the "save" button in a word processor, commits should be made 
early and often.  
A small commit with a precise comment, is often
more useful that a massive commit with a detailed changelist buried in bullets.
If additional comments are necessary, these may be included following a 
blank line as shown below.

```
git commit -m 'Make Brief comments to explain your Feature

Additional detail below a blank line, like this.
-An explanation of why a change was made can help here
-Bullets are useful for a large commit
-But usually, it is better to just commit more often
-With a short one-line description of each little change'
```
____

### 11. Push commits to GitHub
In order to get any accumulated commits in GitHub on your remote branch, 
push the changes. Make sure that you have the correct development branch 
checked out, first.

```
# make sure the right development branch is checked out. There will be "*" next to the checked out branch
git branch -vv

# if the correct development branch is not checked out, change it
git checkout <branch-name>

# push commits to the development branch
git push origin <branch-name>
```
____

### 12. Submit a pull request
After you git push all commits to your branch and you believe 
you are ready to post code to the main codebase, open GitHub 
and issue a pull request. (It will probably be a highlighted 
button at the top of the page -- "New Pull Request"). Of course, 
please test, document, style check 
(see ['Black'](https://pypi.org/project/black/)), 
and prepare (as-needed) demonstration notebooks, etc. 
More information about pull requests can be found 
[here](https://help.github.com/en/github/collaborating-with-issues-and-pull-requests/about-pull-requests). 
From that page:
> "Pull requests let you tell others about changes you've pushed 
to a branch in a repository on GitHub. Once a pull request is 
opened, you can discuss and review the potential changes with 
collaborators and add follow-up commits before your changes 
are merged into the base branch."

A pull request will allow someone else to look at your code with 
you to make sure that it is ready to share with the world. 
Most of the time, someone who was involved with preparing 
the code can be the reviewer; for major changes, it should 
be someone outside the core development team.
____

### 13. Useful references
- If some of the terminology used above is confusing - we agree. 
This [GitHub glossary](https://help.github.com/en/github/getting-started-with-github/github-glossary#checkout) 
is a useful reference.
-  For a better guide on how to use GitHub to collaborate on 
NOAA development projects, please see the 
[documentation associated with the DMOD repo](https://github.com/NOAA-OWP/DMOD/blob/master/doc/GIT_USAGE.md#contributing-tldr) 
