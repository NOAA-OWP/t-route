## Notes:
Cython doesn't understand namespace packages, or if it does, I don't know how to make it use them.  To work around that, in the namespace directory, add `__init__.pxd`.  Python will ignore this and can treat the directory as a namespace package, and Cython will see it as a "normal" package and allow relative imports to work correctly.


Framework and Routing are too dependent to really use two differnt pacakges, even under the same namespace.

It would make sense to refactor the framework so that routing is a subpackage.  Then perhaps implement a subpackage which allows plugging in different network types and routing types?  Some mix of namespace package (for the framework) with plugin like utility.

https://packaging.python.org/guides/creating-and-discovering-plugins/

Should probably move build_tests.py under this routing package...
