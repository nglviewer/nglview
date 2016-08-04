html:
	python  ./scripts/fix_readme_link.py ../nglview/README.md
	cp doc/ngl*.{js,html} dev/
	cp doc/ngl*.{js,html}  latest/
	(cd doc && make html)

clean:
	git clean -fxd

push:
	git push upstream gh-pages

pull:
	git pull upstream gh-pages
