html:
	(cd doc && make html)

clean:
	git clean -fxd

push:
	git push upstream gh-pages
