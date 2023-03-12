import wezel
from wezel.widgets import TableDisplay
from dbdicom.wrappers import sklearn, skimage, scipy, dipy
import webbrowser


def about_menu(parent):

    parent.action(AboutWezel, text='Wezel', icon=wezel.icons.question_mark)
    parent.action(AboutIBEAt, text='iBEAt')  
    parent.action(AboutBEAtDKD, text='BEAt-DKD')


class AboutWezel(wezel.gui.Action):
    def run(self, app):
        webbrowser.open("https://github.com/QIB-Sheffield/wezel")

class AboutIBEAt(wezel.gui.Action):
    def run(self, app):
        webbrowser.open("https://bmcnephrol.biomedcentral.com/articles/10.1186/s12882-020-01901-x")

class AboutBEAtDKD(wezel.gui.Action):
    def run(self, app):
        webbrowser.open("https://www.beat-dkd.eu/")



def segment_menu(parent): 

    parent.action(WholeKidneyMask, text="Presegment whole kidneys")
    parent.action(RenalSinusFat, text="Renal sinus fat")



class WholeKidneyMask(wezel.gui.Action): 

    def run(self, app):

        fat_desc = 'T1w_abdomen_dixon_cor_bh_fat_post_contrast' 
        out_desc = 'T1w_abdomen_dixon_cor_bh_out_phase_post_contrast'

        app.status.message('Finding correct sequences in the list..')

        fat = app.database().series(SeriesDescription=fat_desc)
        out = app.database().series(SeriesDescription=out_desc)
        msg = 'OK'
        if fat==[] or out==[]:
            msg = 'Not all sequences have been found - please check the list'
        elif len(fat)>1 or len(out)>1:
            msg = 'Some sequences have been repeated - please select the appropriate for analysis'
        if msg == 'OK':
            features = fat + out
        else:
            all = app.database().series()
            cancel, features = app.dialog.input(
                {'label': fat_desc, 'type':'select record', 'options': all, 'default': fat},
                {'label': out_desc, 'type':'select record', 'options': all, 'default': out},
                title = msg)
            if cancel:
                return
            
        app.status.message('Preparing for kidney preprocessing..')

        clusters = sklearn.sequential_kmeans(features, n_clusters=2, multiple_series=True)
        app.display(clusters)
        app.refresh()



class RenalSinusFat(wezel.gui.Action):

    def run(self, app):
 
        out_desc = 'T1w_abdomen_dixon_cor_bh_fat_post_contrast'

        app.status.message('Finding correct sequences in the list..')

        out = app.database().series(SeriesDescription=out_desc)
        lk = app.database().series(SeriesDescription='LK')
        rk = app.database().series(SeriesDescription='RK')
        kidneys = lk+rk

        if out==[] or lk==[] or rk==[] or len(out)>1 or len(lk)>1 or len(rk)>1:
            all = app.database().series()
            cancel, selected = app.dialog.input(
                {"label":out_desc, "type":"select record", "options": all, 'default': out},
                {"label":"Whole-kidney regions", "type":"select records", "options": all, 'default': kidneys},
                title = "Please select input for renal sinus fat processing")
            if cancel:
                return
            out = selected[0]
            kidneys = selected[1]
        else:
            out = out[0]

        cleanup = True

        sf_series = []
        fat_image_masked, fat_mask = dipy.median_otsu(out, median_radius=1, numpass=1)
        for kidney in kidneys:
            # Pipeline calculation
            kidney_hull = skimage.convex_hull_image_3d(kidney)
            # sinus = scipy.image_calculator(kidney_hull, kidney, 'series 1 - series 2', integer=True)
            # sinus_fat = scipy.image_calculator(fat_mask, sinus, 'series 1 * series 2', integer=True)
            sinus_fat = scipy.image_calculator(fat_mask, kidney_hull, 'series 1 * series 2', integer=True)
            sinus_fat_largest = scipy.extract_largest_cluster_3d(sinus_fat)
            sinus_fat_largest.SeriesDescription = kidney.instance().SeriesDescription + 'SF'
            # Append and display
            sf_series.append(sinus_fat_largest)
            app.display(sinus_fat_largest, view='Surface')
            # Remove intermediate results
            if cleanup:
                kidney_hull.remove()
                #sinus.remove()
                sinus_fat.remove()
        fat_image_masked.remove()
        if cleanup:
            fat_mask.remove()
        #   Collect features & display
        df = skimage.volume_features(sf_series)
        app.addWidget(TableDisplay(df), 'ROI statistics')
        app.refresh()

