# Guidance on how to contribute

> All contributions to this project will be released to the public domain. By submitting a pull request or filing a bug, issue, or feature request, you are agreeing to comply with this waiver of copyright interest. Details can be found in our [TERMS](TERMS.md) and [LICENSE](LICENSE).


There are two primary ways to help:
 - Using the issue tracker, and
 - Changing the code-base.


## Using the issue tracker

Use the issue tracker to suggest feature requests, report bugs, and ask questions. This is also a great way to connect with the developers of the project as well
as others who are interested in this solution.

Use the issue tracker to find ways to contribute. Find a bug or a feature, mention in the issue that you will take on that effort, then follow the _Changing the code-base_ guidance below.


## Changing the code-base

Generally speaking, you should fork this repository, make a new branch, make changes in your own fork/branch, and then submit a pull request to incorporate changes into the main codebase. All new code *should* have associated unit tests that validate implemented features and the presence or lack of defects. In almost all cases, the submitted code should have some kind of demonstration notebook to go with it so that we can help others experiment with our project.

Additionally, the code should follow any stylistic and architectural guidelines prescribed by the project. In the absence of such guidelines, mimic the styles and patterns in the existing code-base.

## A step-by-step guide to contributing via GitHub

1. On GitHub Create a fork and clone it to your local directory. A step-by-step guide illustrating how to fork and clone a repository can be found [here](https://help.github.com/en/github/getting-started-with-github/fork-a-repo). Once a new fork is created, use this command to clone it locally:

```
git clone github.com/<githubusername>/t-route.git
```

2. Set up the Git configuration parameters to specify your user name and email address: 

```
git config --global user.name "Your GitHub user name"
git config --global user.email "Your email address"
```

3. Set up a remote connection to the newly forked repo. Call the remote "origin" with branch "master"

```
git remote add origin master https://github.com/<user.name>/t-route.git
```

4. If you wish to occasionally fetch and merge updates from the main t-route codebase, set the upstream remote accordingly. Use the `fetch` and `pull` command to occasionally retrieve updates from the main codebase and incorporate them into your local copy. 

```
git remote add upstream https://github.com/NOAA-OWP/t-route.git
git fetch upstream
git merge upstream
```

5. Create your own branch to make targeted contributions. Branches are parallel versions of the repository that allow you to make changes and tinker as much as you'd like without disrupting the master codebase. It is reccomened to create a new branch for each bite-sized contribution you make. Here is how to create a new branch:

```
git branch <new-branch-name>
```

6. Checkout your newly created branch. The `checkout` command helps you navigate between branches. It updates files in your working directory to match the version stored in that branch. 

```
git checkout <branch-name>
```

7. In your own branch, explore the code, make changes, and tinker at will. Try out one of the notebooks or execute one of the `src/python_framework` or `src/python_routing` files -- most have a test script built in.  -- **be creative!**

8. Once you have made the changes you want to make, add changed files to the staging area:

```
git add <list of your changed files>
```

7. Commit your changes to the checked out branch. Commits should be accompanied by a concise and brief comment. Like the "save" button in a wordprocessor, commits should be made early and often.  

```
git commit -m ‘Make Brief comments to explain your Feature
```

8. In order to get any accumulated commits in GitHub on your branch, push the changes:

```
git push
```

9. After you git push all commits to your branch and you believe you are ready to post code to the main codebase, open GitHub and issue a pull request. (It will probably be a highlighted button at the top of the page -- "New Pull Request"). Of course, please test, document, style check (we use ['Black'](https://pypi.org/project/black/), and prepare (as-needed) demonstration notebooks, etc. More information about pull requests can be found [here](https://help.github.com/en/github/collaborating-with-issues-and-pull-requests/about-pull-requests). 
From that page:
> "Pull requests let you tell others about changes you've pushed to a branch in a repository on GitHub. Once a pull request is opened, you can discuss and review the potential changes with collaborators and add follow-up commits before your changes are merged into the base branch."
A pull request will allow someone else to look at your code with you to make sure that it is ready to share with the world. Most of the time, someone who was involved with preparing the code can be the reviewer; for 
major changes, it should be someone outside the core development team.

10. If some of the terminology used above is confusing - we agree. This [GitHub glossary](https://help.github.com/en/github/getting-started-with-github/github-glossary#checkout) is a useful reference.
