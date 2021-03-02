#ifdef GUI 
#include <QtGui/QApplication>
#include <QDialog>
#include <qapplication.h>
#include <QSplashScreen>
#include <Qt>
#include <QPixmap>
#include <QObject>
#include <QThread>
#include "mainwindow.h"
#endif

#include "common.h"
#include "inrichmain.h"
#include "aligatormain.h"
#include "inrichmultimain.h"

#include <iostream>
#include <cstdlib>
#include <stdio.h>


int main(int argc, char *argv[])
{
    #ifdef GUI 
        QApplication app(argc, argv);

        MainWindow *w = new MainWindow(QString::fromStdString(inrich_ver));
        w->setWindowTitle( QString::fromStdString(inrich_ver) + ": Interval-based Enrichment Analysis for Genome Wide Association Studies");
        w->setWindowIcon(QIcon(":/images/inrich_log.jpg"));
        w->show();

        return app.exec();
    #endif
    
    #ifndef GUI 
    if(argc>1){
        InputData globals;
        globals.set_data(argc, argv);
        std::string msg = globals.check_input_data();
        if(msg.compare("")==0) {
            if(globals.multi_set) {
                InRichMultiMain inrich(globals);                
                if(globals.exam_bp_dist) {
                    std::cerr << "Physical Clustering Test cannot be run on the Multi Set data\n\n";
                    return 1;
                }
                if(globals.convert_mode) {
                    std::cerr << "SNP2Interval Conversion cannot be run with Multi Set Analysis\n\n";
                    return 1;
                }
                else {
                    inrich.inrich_multi_test();
                }
            }
            else if(!globals.aligator) {
                InRichMain inrich(globals);

                if(globals.exam_bp_dist) {
                    inrich.clustering_test();
                    return 0;
                }

                inrich.inrich_test();
            }
	    else {
                AligatorMain aligator;
                aligator.set_globals(globals);
                aligator.aligator_test();
	    }
        }
        else {
            std::cout << msg;
        }
	return 0;
    }

    #endif
}
