import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib_venn import venn2, venn2_circles

import brewer2mpl
bmap = brewer2mpl.get_map('Set1', 'qualitative', 9)
mypal = bmap.mpl_colors

# plot GM12878 prediction and experiment venn diagram
plt.figure(figsize=(4.2,3.6))

v = venn2(subsets=(12,340,120),set_labels=None)
plt.setp(v.get_patch_by_id('10'), 'color', 'white')
plt.setp(v.get_patch_by_id('11'), 'color', 'white')
plt.setp(v.get_patch_by_id('01'), 'color', 'white')

value_10 = plt.getp(v.get_label_by_id('10'),'text')
plt.setp(v.get_label_by_id('10'),  'text','')

c = venn2_circles(subsets=(12,340,120))
plt.setp(c,
    'linestyle', 'dashed')
plt.setp(c[0],
    'color', mypal[0],
    'fill',False)

plt.title('GM12878')
plt.annotate(
    'Experiment', 
    xy=v.get_label_by_id('10').get_position(), 
    xytext= (-30,35),
    color=mypal[0],
    ha='center', 
    textcoords='offset points' )
plt.annotate(
    'Prediction', 
    xy=v.get_label_by_id('01').get_position(),
    xytext=(65/3,90),
    color="black",
    ha='center', 
    textcoords='offset points' )
plt.annotate(
    value_10,
    xy=v.get_label_by_id('10').get_position() + np.array([0.01, 0]), 
    xytext=(-40,-40),
    textcoords='offset points',arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0.5',color='gray'))
plt.savefig('/tmp/GM.Pred.Exp.Venn.Sig.eps')

# ------------------------------------------------
# plot K562 prediction and experiment venn diagram
plt.figure(figsize=(4.2,3.6))

v = venn2(subsets=(3,346,80),set_labels=None)
plt.setp(v.get_patch_by_id('10'), 'color', 'white')
plt.setp(v.get_patch_by_id('11'), 'color', 'white')
plt.setp(v.get_patch_by_id('01'), 'color', 'white')

value_10 = plt.getp(v.get_label_by_id('10'),'text')
plt.setp(v.get_label_by_id('10'),  'text','')

c = venn2_circles(subsets=(3,346,80))
plt.setp(c,
    'linestyle', 'dashed')
plt.setp(c[0],
    'color', mypal[1],
    'fill',False)

plt.title('K562')
plt.annotate(
    'Experiment', 
    xy=v.get_label_by_id('10').get_position(), 
    xytext= (-30,35),
    color=mypal[1],
    ha='center', 
    textcoords='offset points' )
plt.annotate(
    'Prediction', 
    xy=v.get_label_by_id('01').get_position(),
    xytext=(65/3,90),
    color="black",
    ha='center', 
    textcoords='offset points' )
plt.annotate(
    value_10,
    xy=v.get_label_by_id('10').get_position() + np.array([0.01, 0]),
    xytext=(-40,-40),
    textcoords='offset points',arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0.5',color='gray'))
plt.savefig('/tmp/K.Pred.Exp.Venn.Sig.eps')

# ------------------------------------------------
# plot CHIA-PET prediction and experiment venn diagram
plt.figure(figsize=(4.2,3.6))

v = venn2(subsets=(3,0,10),set_labels=None)
plt.setp(v.get_patch_by_id('10'), 'color', 'white')
plt.setp(v.get_patch_by_id('11'), 'color', 'white')
plt.setp(v.get_patch_by_id('01'), 'color', 'white')

value_10 = plt.getp(v.get_label_by_id('10'),'text')
plt.setp(v.get_label_by_id('10'),  'text','')
value_01 = plt.getp(v.get_label_by_id('11'),'text')
plt.setp(v.get_label_by_id('01'),  'text','')

c = venn2_circles(subsets=(3,0,10))
plt.setp(c,
    'linestyle', 'dashed')
plt.setp(c[0],
    'color', mypal[1],
    'fill',False)

plt.title('K562')
plt.annotate(
    'CHIA-PET', 
    xy=v.get_label_by_id('10').get_position(), 
    xytext= (-30,35),
    color=mypal[1],
    ha='center', 
    textcoords='offset points' )
plt.annotate(
    'Prediction', 
    xy=c[1].center,
    xytext=(85/3,35),
    color="black",
    ha='center', 
    textcoords='offset points' )
plt.annotate(
    value_10,
    xy=v.get_label_by_id('10').get_position() + np.array([0.01, 0]),
    xytext=(-40,-40),
    textcoords='offset points',arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0.5',color='gray'))
plt.savefig('/tmp/K.Pred.Exp.Venn.CHIA-PET.eps')
