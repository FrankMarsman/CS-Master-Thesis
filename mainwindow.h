#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <sim2d.h>
#include <ctime>
#include <QVector>
#include <QPen>

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
  Q_OBJECT

public:
  Sim2D * sim;
  bool isRunning; // true if NextStep is called continuously

  double simFPS; // how many times NextStep per sec

  uint prevRunStep, curRunStep; // to comute net FPS
  time_t simStepTimer;  // to comute net FPS
  double netFPS; // fps of everything; nextstep, update screen etc

  QVector <double> E_tot, E_pot, E_grav, E_kin;
  QVector <double> tVec;

  // for experiment:
  double SIM_MAX_RUNTIME;
  QVector <double> expVecX;
  QVector < QVector <double> > expVecYs; // all lines
  QVector <QString> yVecLabels; // for legend
  QVector <QPen> yVecPens; // pen for each line

  explicit MainWindow(QWidget *parent = 0);
  ~MainWindow( );

public slots:
  void Run( );
  void UpdateInfoText( ); // updates information on ui
  void UpdateMeshImage( ); // prints mesh to ui
  void UpdateSelectedVertices( );
  void PerformInitDeform( );

  void UpdateStepPlot( );

  void RunExperiment( );

private slots:
  void on_nextButton_clicked( );
  void on_noiseButton_clicked( );

  void on_initButton_clicked();

  void on_MainTab_currentChanged(int index);

  void on_addVButton_clicked();


  void on_lockButton_clicked();

  void on_resetImgButton_clicked();

  void on_vNoiseButton_clicked();

  void on_stepButton_clicked();

  void on_squeezeXButton_clicked();


  void on_squeezeYButton_clicked();

  void on_expButton_clicked();

  void on_copyExpBut_clicked();

  void on_copyMeshImgButton_clicked();

private:
  Ui::MainWindow *ui;
};

#endif // MAINWINDOW_H
