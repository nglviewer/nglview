# Example

```python
import seaborn as sns

view = nglview.show_mdtraj(traj)
view.clear()
view

view.add_cartoon('protein')
# get color palette from seaborn
sns_color_palette = sns.cubehelix_palette(traj.n_residues).as_hex()

# convert to int
colors = [int(c.replace('#', '0x')) for c in sns_color_palette]
view._set_color_by_residue(colors)
```
