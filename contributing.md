# Guidance on how to contribute

> All contributions to this project will be released to the public domain.
> By submitting a pull request or filing a bug, issue, or
> feature request, you are agreeing to comply with this waiver of copyright interest.
> Details can be found in our [TERMS](TERMS.md) and [LICENSE](LICENSE).


There are two primary ways to help:
 - Using the issue tracker, and
 - Changing the code-base.


## Using the issue tracker

Use the issue tracker to suggest feature requests, report bugs, and ask questions.
This is also a great way to connect with the developers of the project as well
as others who are interested in this solution.

Use the issue tracker to find ways to contribute. Find a bug or a feature, mention in
the issue that you will take on that effort, then follow the _Changing the code-base_
guidance below.


## Changing the code-base

Generally speaking, you should fork this repository, make changes in your
own fork, and then submit a pull request. All new code *should* have associated
unit tests that validate implemented features and the presence or lack of defects.
In almost all cases, the submitted code should have some kind of demonstration notebook
to go with it so that we can help others experiment with our project.

Additionally, the code should follow any stylistic and architectural guidelines
prescribed by the project. In the absence of such guidelines, mimic the styles
and patterns in the existing code-base.

## The process in summary
1. On GitHub Create a Fork
`git clone github.com/<githubusername>/t-route.git`
2. Try out one of the notebooks or execute one of the `src/python_framework` or `src/python_routing files` -- most have a test script built in.
3. Make your changes to the source files -- **be creative!**
4. If you are working on a complex feature that may require multiple days of effort, consider making a branch and committing to that branch.
5. If you are working with someone else to develop a feature, you can share the fork with them and follow these same steps on the fork until
you are jointly ready to submit to the main repository.
6. Once you have made the changes you want to make, add ...
```
git add <list of your changed files>
```
7. ... and commit them. As suggested in the example, a small commit with a precise comment, is often
more useful that a massive commit with a detailed changelist buried in bullets.
```
git commit -m ‘Make Brief comments to explain your Feature

If you need additional lines for more detail, add them below a double space, like this.
-Bullets are useful for a large commit
-But usually, it’s better to just commit more often
-With a short one-line description of each little change’
```
8. In order to get any accumulated commits in GitHub on your fork, push the changes using `git push`
9. After you git push all commits to your fork and you believe you are ready to post code to the main 
repository, open GitHub and issue a pull request. (It will probably be a highlighted button at the top of the 
page -- "New Pull Request").
Of course, please test, document, style check (we use ['Black'][https://pypi.org/project/black/]), and prepare (as-needed) demonstration notebooks, etc. 
More information about pull requests can be found [here][https://help.github.com/en/github/collaborating-with-issues-and-pull-requests/about-pull-requests]. 
From that page:
> "Pull requests let you tell others about changes you've pushed to a branch in a repository on GitHub. 
Once a pull request is opened, you can discuss and review the potential changes with collaborators and 
add follow-up commits before your
changes are merged into the base branch."
A pull request will allow someone else to look at your code with you to make sure that it is ready to share 
with the world. Most of the time, someone who was involved with preparing the code can be the reviewer; for 
major changes, it should be someone outside the core development team.
9. **IMPORTANT** After you have issued a pull request the master upstream repository (NOAA-OWP) will have 
been updated with your new code.
It is important to make sure your fork is kept up-to-date with these new changes with the following commands. 
The first command
adds the OWP repository as an upstream remote that can be 'fetched' from to get the updates.
```
git remote add upstream https://github.com/NOAA-OWP/t-route.git
git fetch upstream && git merge upstream/master

```
