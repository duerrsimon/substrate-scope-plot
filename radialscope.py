from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import Draw
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib import rcParams
import os
from svgutils import transform,compose

class strSVG(compose.Element):
    """SVG from string.

    Parameters
    ----------
    fname : str
       full path to the file
    """

    def __init__(self, svg):
        obj = transform.fromstring(svg)
        self.root = obj.getroot().root
        
from IPython.display import SVG # /!\ note the 'SVG' function also in svgutils.compose
import numpy as np
from rdkit import Geometry


def draw_with_indeces(settings):
    """
    Drawing function that displays the input smiles string with all atom indeces
    """
    m = Chem.MolFromSmiles(settings['SMILESSTRING'])
    dm = Draw.PrepareMolForDrawing(m)
    d2d = Draw.MolDraw2DSVG(350,350)
    opts = d2d.drawOptions()
    for i in range(m.GetNumAtoms()):
        opts.atomLabels[i] = m.GetAtomWithIdx(i).GetSymbol()+str(i)
    d2d.DrawMolecule(dm)
    d2d.FinishDrawing()
    return d2d.GetDrawingText()

class RadialScope(object):
    """
        Radial Scope Plot using Matplotlib and RDKIT

        Written by Simon Duerr with code from Greg Landrum for the automatic positioning (@rdkit)

        License: MIT 

        This class handles all the heavy work and constructs a colorbar svg and a pie chart svg for each passed radial scope dictionary and then assembles the svg. 

        The size of the matplotlib plot is always the same, which is why the hard coded center used for positioning on the atoms is likely fine.
    """

    def __init__(self, settings_dict,*args):
        """
        Calls the main constructur. 


        Parameters
        ----------
        settings_dicte : dict
               main settings
           args: dict
               dictionaries containing the information for the scope plots
        """
        self.settings=settings_dict
        self.plots=[]
        for plot in args:
            self.plots.append(plot)
        if len(args)<1:
            raise Exception('minimum 1 radial scope dictionary must be given')
        self.main()


    def draw_smiles(self):
        """
            draws the smiles string for the plotting, needs to return the drawing objects also because the atom indeces are extracted.
        """
        m = Chem.MolFromSmiles(self.settings['SMILESSTRING'])
        dm = Draw.PrepareMolForDrawing(m)
        d2d = Draw.MolDraw2DSVG(350,350)
        if self.settings['use_bw_atom_theme']:
            d2d.drawOptions().useBWAtomPalette()
        d2d.DrawMolecule(dm)
        d2d.FinishDrawing()
        return d2d.GetDrawingText(), d2d, dm

    def replace_label_with_smiles(self,svg_file='', smiles='C=C', search_index='~0'):
        """
        draws a small organic rest, looks for the ~index comment in the complete svg, finds the glyphs and replaces them with the organic subsituent. 
        Note that this is hacky and the position of the organic subsituent likely needs to be fixed in a vector software such as Inkscape. 

        Parameters
        ----------
        svg_file : str
           the text of the svg file which contains the comments that are replaced
        search_index: str
            the name of the comment which is indicated by a prepended tilde by the user
        smiles: str
            the smiles string used for the replacement
        """
        m = Chem.MolFromSmiles(smiles)
        dm = Draw.PrepareMolForDrawing(m)
        d2d = Draw.MolDraw2DSVG(100,100)
        d2d.drawOptions().padding=0
        d2d.drawOptions().clearBackground=False
        d2d.drawOptions().useBWAtomPalette()
        d2d.DrawMolecule(dm)
        d2d.FinishDrawing()
        # scale smiles molecule and remove clutter
        group1 = d2d.GetDrawingText()
        replace_str=group1[group1.find('<!-- END OF HEADER -->')+len("<!-- END OF HEADER -->")+1:-8]
        replace_str='<g transform="translate(-300,-300)scale(6)">'+replace_str+"</g>"
        # find the index in the pie chart that needs to be replaced, we will geplace the two glyphs with the svg text from rdkit
        index_of_comment=svg_file.find(str(search_index))
        index_of_defsend=svg_file[index_of_comment:].find('</defs>')
        start=svg_file[index_of_comment+index_of_defsend:].find('<use')
        end=svg_file[index_of_comment+index_of_defsend:].find('</g>')
        item_to_replace=svg_file[index_of_comment+index_of_defsend+start:index_of_comment+index_of_defsend+end]
        return svg_file.replace(item_to_replace, replace_str)
        
    def plot_figure_and_colorbar(self,radial_scope_setup, vals):
        """
        plots one pie chart and the corresponding colorbar

        Parameters
        ----------
        radial_scope_setup : dict
           contains information about this pie chart like labels, colormap etc.
        vals :list
            contains 5 list with prepended input for the empty wedge in the pie chart and the processed labels.

        Returns
        -------
        fig_svg: str
            the svg for the pie plot
        colorbars
            the svg for the colorbar of the inner and outer circle

        """
        
        fig, ax = plt.subplots(1,figsize=(10,10))

        size = 0.5 # size of inner plot
        alpha = 0
        which_wedge = 0 # first wedge is always transparent
        _=ax.set(aspect="equal")
        circle1 = plt.Circle((0, 0), 0.15, color='w', ls='-', ec='k', lw=1.4, zorder=99, gid='circle_anchor')
        label = ax.annotate(radial_scope_setup['rest_label'], xy=(0, 0), fontsize=30, ha="center", va='center', zorder=100, gid='circle_content')

        cmap_inner = plt.get_cmap(radial_scope_setup['CMAPINNER'])
        cmap_outer = plt.get_cmap(radial_scope_setup['CMAPOUTER'])
        norm_outer = plt.Normalize(np.min(vals[2]), np.max([vals[2]]))
        outer_colors = cmap_outer(norm_outer(vals[2]))
        norm_inner = plt.Normalize(np.min(vals[1]), np.max([vals[1]]))
        inner_colors = cmap_inner(norm_inner(vals[1]))

        labels_circle=ax.pie(vals[0], startangle=radial_scope_setup['startangle'], radius=0.7, colors=['w']*len(vals[0]), labels=vals[5], labeldistance=1.1,
               wedgeprops=dict(width=size, edgecolor='k',linewidth= 1.4), textprops=dict(fontsize='large', weight="semibold",va='center', ha='center') )


        outer_circle=ax.pie(vals[0], radius=0.7,startangle=radial_scope_setup['startangle'], colors=outer_colors, labels=vals[3], labeldistance=0.5,
               wedgeprops=dict(width=size, edgecolor='k',linewidth= 1.4),textprops=dict(fontsize='large',weight="semibold",va='center', ha='center'))

        inner_circle=ax.pie(vals[0], startangle=radial_scope_setup['startangle'], radius=1-size, colors=inner_colors, labels=vals[4], labeldistance=1.2,
               wedgeprops=dict(width=size, edgecolor='k',linewidth= 1.4,), textprops=dict(fontsize='large', weight="semibold",va='center', ha='center') )



        inner_circle[0][which_wedge].set_alpha(alpha)
        outer_circle[0][which_wedge].set_alpha(alpha)
        labels_circle[0][which_wedge].set_alpha(alpha)

        for i in range(len(inner_circle[1])):
            if vals[2][i]>self.settings['white_cutoff']:
                _=inner_circle[1][i].set_color('white')
            if vals[1][i]>self.settings['white_cutoff']:
                _=outer_circle[1][i].set_color('white')
        _=ax.add_artist(circle1)

        from io import StringIO
        
        sio = StringIO()
        _=fig.savefig(sio, transparent=True, format='SVG')
        fig_svg = sio.getvalue()
        _=plt.close(fig)
        # close first figure as we do not want to display it and start assembling the colorbar
        

        fig, axs = plt.subplots(2, figsize=(3,2))
        sm = ScalarMappable(cmap=cmap_inner, norm=plt.Normalize(0,max(vals[1])))
        _=sm.set_array([])
        cbar = plt.colorbar(sm,cax=axs[0], orientation="horizontal",ticks=[0,50,99])
        _=cbar.ax.set_xticklabels(['0', '50', '100 %']) 
        _=cbar.set_label(radial_scope_setup['INNERLABEL'],weight='bold', fontsize=12)
        sm1 = ScalarMappable(cmap=cmap_outer, norm=plt.Normalize(0,max(vals[1])))
        _=sm1.set_array([])
        cbar1 = plt.colorbar(sm1,cax=axs[1], orientation="horizontal",ticks=[0,50,99])
        _=cbar1.ax.set_xticklabels(['0', '50', '100 %']) 
        _=cbar1.set_label(radial_scope_setup['OUTERLABEL'],weight='bold',fontsize=12)
        figure2=fig.tight_layout()
        sio2 = StringIO()
        _=fig.savefig(sio2, transparent=True, format='SVG')
        colorbars = sio2.getvalue()
        _=plt.close(fig)
        return fig_svg, colorbars

    def main(self):
        settings=self.settings
        SMILESSTRING=settings['SMILESSTRING']
        resulting_plots=[]
        pRList=[]
        mol_svg, d2d, dm=self.draw_smiles()
        replace_index=[]
        for scope_plot in self.plots:
            # for each scope plot, make a vals list containing empty first items for the wedge with alpha=0
            if type(scope_plot)!=dict:
                continue
            
            sizes=[360-scope_plot['coverangle_wedges']]+[scope_plot['coverangle_wedges']/scope_plot['no_wedges']]*scope_plot['no_wedges']
            
            label_inner_circle, label_outer_circle=['']+['']*scope_plot['no_wedges'],['']+['']*scope_plot['no_wedges']
            if (len(scope_plot['value_inner_circle'])!=scope_plot['no_wedges'] or len(scope_plot['value_outer_circle'])!=scope_plot['no_wedges']):
                print('not equal')
            value_inner_circle, value_outer_circle=scope_plot['value_inner_circle'], scope_plot['value_outer_circle']
            rounding_boundary=scope_plot['rounding_boundary']
            value_groups=scope_plot['value_groups']
            
            for i in range(scope_plot['no_wedges']):
                if scope_plot['rounding']:
                    if value_inner_circle[i]>=rounding_boundary:
                        label_inner_circle[i+1]=">"+str(value_inner_circle[i])
                    else: 
                        label_inner_circle[i+1]=str(value_inner_circle[i])

                    if value_outer_circle[i]>=rounding_boundary:
                        label_outer_circle[i+1]=">"+str(value_outer_circle[i])
                    else:
                        label_outer_circle[i+1]=str(value_outer_circle[i])
                else:
                    label_inner_circle[i+1]=str(value_inner_circle[i])
                    label_outer_circle[i+1]=str(value_outer_circle[i])
            j=0
            for i,item in enumerate(value_groups):
                if item[0]=='~':
                    replace_index.append(('~'+str(j),item[1:]))
                    value_groups[i]='~'+str(j)
                    j=j+1

            vals = [sizes, # size of the wedges, the first wedge is transparent and will not be shown 
                    [0]+value_inner_circle, # colormap values for the inner circle, maximum value determines intensity, first is for the transparent wedge and should stay 0
                    [0]+value_outer_circle, # colormap values for the outer circle, maximum value determines intensity, first is for the transparent wedge and should stay 0
                    label_inner_circle, #labels for the inner circle
                    label_outer_circle, #labels for the outer circle    
                    [""]+value_groups, #groups  
                   ]
            resulting_plots.append(self.plot_figure_and_colorbar(scope_plot,vals))

            # get the atom id from the settings and save its position
            rIdx = scope_plot['attach_atom_id']
            pRList.append(d2d.GetDrawCoords(Geometry.Point2D(dm.GetConformer().GetAtomPosition(rIdx))))
        
        # take colorbar from first plot  #ToDo extension to multiple colorbars
        colorbar=compose.Panel(strSVG(resulting_plots[0][1]).scale(0.8).move(-350,400))
        panels=[compose.Panel(strSVG('<svg></svg>'))]*len(resulting_plots)
        for i,plot in enumerate(resulting_plots):
            panels[i]=strSVG(resulting_plots[i][0]).move(-369,-358).scale(1).move(pRList[i].x,pRList[i].y)
        compose.Figure("720", "720", 
            compose.Panel(strSVG(mol_svg).scale(1).move(0,0)),
            colorbar,
            *panels
            ).move(350,350).scale(self.settings['scalefactor']).save("substrate_scope.svg")
        new_svg=SVG('substrate_scope.svg')._data
        for item in replace_index:
            new_svg=self.replace_label_with_smiles(svg_file=new_svg, smiles=item[1],search_index= item[0] )

        if settings['use_bold_font']:
            new_svg.replace('font-weight:normal', 'font-weight:bold')
        f = open("substrate_scope_replaced.svg","w") 
        f.write(new_svg)
        f.close()
        print('File written to:',os.getcwd()+'/substrate_scope_replaced.svg')
        
        
        

        
        