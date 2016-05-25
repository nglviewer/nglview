html:
	(cd doc && make html)

clean:
	git clean -fxd

push:
	git commit -m 'update doc'
	git push upstream gh-pages

pull:
	git pull upstream gh-pages
