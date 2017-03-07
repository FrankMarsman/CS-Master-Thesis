#include "mainwindow.h"
#include "ui_mainwindow.h"

#include <iostream>
#include <Dense>
#include <Eigen>
#include <QDebug>
#include <QImage>
#include <QTimer>
#include <ctime>
#include <QString>
#include <cmath>
#include <QVector>

#include <qcustomplot.h>

#include <QMessageBox>
#include <QApplication>

using namespace Eigen;
using namespace std;

MainWindow::MainWindow(QWidget *parent) :
  QMainWindow(parent),
  ui(new Ui::MainWindow) {
  ui->setupUi(this);

  uint N = 15;
  this->sim = new Sim2D(N * N);

  // lock bottom vertices
  sim->SetSelectedVertices(0, 1, 0, (2.0 / N));
  sim->AddLockedVertices( );
  sim->selectedVertices.clear( );

  // make sure ui updates when settings change:
  QObject::connect(ui->aaBox, SIGNAL(toggled(bool)), SLOT(UpdateMeshImage( )));
  QObject::connect(ui->drawEdgeBox, SIGNAL(toggled(bool)), SLOT(UpdateMeshImage( )));
  QObject::connect(ui->sX1Slider, SIGNAL(actionTriggered(int)), SLOT(UpdateSelectedVertices( )));
  QObject::connect(ui->sX2Slider, SIGNAL(actionTriggered(int)), SLOT(UpdateSelectedVertices( )));
  QObject::connect(ui->sY1Slider, SIGNAL(actionTriggered(int)), SLOT(UpdateSelectedVertices( )));
  QObject::connect(ui->sY2Slider, SIGNAL(actionTriggered(int)), SLOT(UpdateSelectedVertices( )));

  prevRunStep = 0;
  simStepTimer = clock( );





  ui->mBox->setValue(N);

  isRunning = true; // set to opposite what you want

  on_nextButton_clicked( );
  on_resetImgButton_clicked( );

  UpdateSelectedVertices( );
  UpdateMeshImage( );

  QTimer::singleShot(100, this, SLOT(Run( )));

  // maximize window:
  this->setWindowState(Qt::WindowMaximized);
} // MainWindow

MainWindow::~MainWindow( ) {
  delete ui;
} // ~MainWindow

void MainWindow::Run( ) {

  if (isRunning && this->sim->h != ui->hValSpinBox->value( )) {
    this->sim->h = ui->hValSpinBox->value( );

    if (this->sim->simStep != 0) // because if so sim will init lMatrix anyway
      this->sim->InitLMatrix( );
  } // if

  this->sim->iterationsPerStep = ui->nStepsSpinBox->value( );

  bool useLLT = ui->lltBox->isChecked( );

  if (isRunning) {
    bool enableGrav = ui->gravBox->isChecked( );
    double g = (enableGrav)? 0.01 : 0;
    sim->SetGravity(g);

    double C = ui->sfBox->value( ); // spring force constant
    if (sim->springForceConstant != C)
      sim->SetSpringForceConstant(C);

    sim->doParallelPi = ui->parBox->isChecked( );
    sim->doParallelRmat = ui->parRBox->isChecked( );

    uint numSteps = ui->frameStepBox->value( );

    time_t t = clock( );

    for (uint i = 0; i < numSteps; i++)
      this->sim->NextStep(useLLT);

    // compute stats:
    double elapsed = double(clock( ) - t) / CLOCKS_PER_SEC;
    double curFPS = double(numSteps) / elapsed;
    this->simFPS = (simFPS * 0.9) + 0.1 * curFPS;

    if (isinf(simFPS))
      simFPS = 0;

    if (double(clock( ) - simStepTimer) >= double(CLOCKS_PER_SEC)) {
      elapsed = double(clock( ) - simStepTimer) / CLOCKS_PER_SEC;
      this->netFPS = double(curRunStep - prevRunStep) / elapsed;
      prevRunStep = curRunStep;
      simStepTimer = clock( );
    } // if

    UpdateMeshImage( );
    UpdateInfoText( );

    curRunStep++;
  } // if
  QTimer::singleShot(1, this, SLOT(Run( )));
} // Run


void MainWindow::UpdateInfoText( ) {
  ui->fpsLabel->setText(QString::number(simFPS, 'f', 0));
  ui->tLabel->setText(QString::number(sim->simStep, 'e', 2));
  ui->netFpsLabel->setText(QString::number(netFPS, 'f', 0));
  ui->simInfoLabel->setText(sim->GetInfoString( ));

  // TEMP!
  QCustomPlot * customPlot = ui->ePlot;

  this->tVec.push_back(sim->simStep);
  this->E_grav.push_back(sim->GetGravPotEnergy( ));
  this->E_pot.push_back(sim->GetSpringPotentialEnergy( ));
  this->E_kin.push_back(sim->GetKineticEnergy( ));
  this->E_tot.push_back(E_grav.back( ) + E_kin.back( ) + E_pot.back( ));

  // create graph and assign data to it:
  customPlot->clearGraphs( );

  customPlot->addGraph( );
  customPlot->addGraph( );
  customPlot->addGraph( );
  customPlot->addGraph( );

  //qDebug( ) << "number of grapgs:" << customPlot->graphCount( );

  customPlot->graph(0)->setData(tVec, E_tot);
  customPlot->graph(1)->setData(tVec, E_pot);
  customPlot->graph(2)->setData(tVec, E_kin);
  customPlot->graph(3)->setData(tVec, E_grav);

  customPlot->graph(0)->setPen(QPen(Qt::blue));
  customPlot->graph(1)->setPen(QPen(Qt::red));
  customPlot->graph(2)->setPen(QPen(Qt::green));
  customPlot->graph(3)->setPen(QPen(Qt::black));

  // give the axes some labels:
  customPlot->xAxis->setLabel("x");
  customPlot->yAxis->setLabel("y");
  // set axes ranges, so we see all data:
  customPlot->xAxis->setRange(0, tVec.back( ));
  customPlot->yAxis->setRange(0, E_tot.back( ) * 2.0);

  customPlot->replot( );
} // UpdateInfoText

void MainWindow::UpdateMeshImage( ) {

  bool useAA = ui->aaBox->isChecked( );
  bool drawEdges = ui->drawEdgeBox->isChecked( );
  bool drawSelectedVertices = (ui->MainTab->currentIndex( ) == 1);
  bool drawLockedVertices = (ui->MainTab->currentIndex( ) == 0);

  uint DRAW_MODE = ui->drawModeBox->currentIndex( );

  //qDebug( ) << "UpdateMeshImage, DRAW_MODE =" << DRAW_MODE;

  QImage img = sim->ToQImage(500, useAA, drawEdges, drawSelectedVertices, drawLockedVertices, DRAW_MODE);
  ui->imgLabel->setPixmap(QPixmap::fromImage(img));
} // UpdateMeshImage

void MainWindow::UpdateSelectedVertices( ) {

  double x1 = double(ui->sX1Slider->value( )) / 100;
  double x2 = double(ui->sX2Slider->value( )) / 100;
  double y2 = double(ui->sY1Slider->value( )) / 100;
  double y1 = double(ui->sY2Slider->value( )) / 100;

  this->sim->SetSelectedVertices(x1, x2, 1.0 - y1, 1.0 - y2);
  UpdateMeshImage( );
} // UpdateSelectedVertices

// adds velocity to selected vertices
void MainWindow::on_addVButton_clicked( ) {
  if (sim->selectedVertices.size( ) == 0) { // no vertices to add noise to
    QMessageBox::information(
        this,
        tr("No!"),
        tr("No vertices selected! So can't add velocity.") );
    return;
  } // if

  double vx = ui->vxBox->value( );
  double vy = ui->vyBox->value( );

  this->sim->AddVelocity(vx, vy);
  UpdateMeshImage( );
} // on_addVButton_clicked

// compute next state of sim and show on screen
void MainWindow::on_nextButton_clicked( ) {
  if (isRunning) {
    isRunning = false;
    ui->nextButton->setText("Start running");
  } // if
  else {
    isRunning = true;
    ui->nextButton->setText("Stop running");
  } // else
} // on_nextButton_clicked


void MainWindow::on_noiseButton_clicked( ) {
  if (sim->selectedVertices.size( ) == 0) { // no vertices to add noise to
    QMessageBox::information(
        this,
        tr("No!"),
        tr("No vertices selected! So can't add noise.") );
    return;
  } // if
  this->sim->AddNoise(0.1);
  UpdateMeshImage( );
} // on_noiseButton_clicked

// initializes this->sim with user given number of vertices
void MainWindow::on_initButton_clicked( ) {
  if (isRunning) { // sim is running; dangerous; return
    QMessageBox::information(
        this,
        tr("No!"),
        tr("Can't reset sim because it is still running! Stop first.") );
    return;
  } // if

  uint m = ui->mBox->value( );
  m *= m; // make squared

  delete this->sim;
  this->sim = NULL;
  this->sim = new Sim2D(m);

  on_resetImgButton_clicked( );

  this->tVec.clear( );
  this->E_grav.clear( );
  this->E_pot.clear( );
  this->E_kin.clear( );
  this->E_tot.clear( );

  UpdateSelectedVertices( );
  UpdateMeshImage( );
} // on_initButton_clicked

// user selects different tab
void MainWindow::on_MainTab_currentChanged(int index) {
  qDebug( ) << "Tab changed to" << index;
  UpdateMeshImage( );
} // on_MainTab_currentChanged

// locks positions of currently selected vertices
void MainWindow::on_lockButton_clicked( ) {
  if (sim->selectedVertices.size( ) == 0) {
    QMessageBox::information(
        this,
        tr("No!"),
        tr("No vertices selected to lock!") );
    return;
  } // if

  this->sim->AddLockedVertices( );
} // on_lockButton_clicked


// recenteres view
void MainWindow::on_resetImgButton_clicked( ) {
  this->sim->InitImgParams( );

  ui->imgWBox->setValue(sim->imgViewSize);
  ui->imgXBox->setValue(sim->imgCenterX);
  ui->imgYBox->setValue(sim->imgCenterY);

  UpdateMeshImage( );
} // on_resetImgButton_clicked


// adds random velocity to selected vertices
void MainWindow::on_vNoiseButton_clicked( ) {
  this->sim->AddVNoise( );

  UpdateMeshImage( );
} // on_vNoiseButton_clicked

















