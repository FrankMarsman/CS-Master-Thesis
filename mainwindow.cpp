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

#include <QWinTaskbarButton>
#include <QWinTaskbarProgress>

using namespace Eigen;
using namespace std;

MainWindow::MainWindow(QWidget *parent) :
  QMainWindow(parent),
  ui(new Ui::MainWindow) {
  ui->setupUi(this);

  uint N = 20;
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


  // init pens for plots:
  this->yVecPens.push_back(QPen(Qt::blue));
  this->yVecPens.push_back(QPen(Qt::red));
  this->yVecPens.push_back(QPen(Qt::green));
  this->yVecPens.push_back(QPen(Qt::yellow));
  this->yVecPens.push_back(QPen(Qt::magenta));

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


// runs a single step in simulation
void MainWindow::on_stepButton_clicked( ) {
  if (this->sim->h != ui->hValSpinBox->value( )) {
    this->sim->h = ui->hValSpinBox->value( );

    if (this->sim->simStep != 0) // because if so sim will init lMatrix anyway
      this->sim->InitLMatrix( );
  } // if

  bool useLLT = ui->lltBox->isChecked( );

  this->sim->iterationsPerStep = ui->nStepsSpinBox->value( );
  // gravity:
  bool enableGrav = ui->gravBox->isChecked( );
  double g = (enableGrav)? 0.1 : 0;
  if (sim->gravAcc != g)
    sim->SetGravity(g);

  // spring constant:
  double C = ui->sfBox->value( ); // spring force constant
  if (sim->springForceConstant != C)
    sim->SetSpringForceConstant(C);

  this->sim->NextStep(useLLT);

  UpdateMeshImage( );
  UpdateInfoText( );

  curRunStep++;
} // on_stepButton_clicked

void MainWindow::Run( ) {

  if (isRunning && this->sim->h != ui->hValSpinBox->value( )) {
    this->sim->h = ui->hValSpinBox->value( );

    if (this->sim->simStep != 0) // because if so sim will init lMatrix anyway
      this->sim->InitLMatrix( );
  } // if

  this->sim->iterationsPerStep = ui->nStepsSpinBox->value( );

  bool useLLT = ui->lltBox->isChecked( );

  if (isRunning) {

    // gravity:
    bool enableGrav = ui->gravBox->isChecked( );
    double g = (enableGrav)? 0.1 : 0;
    if (sim->gravAcc != g)
      sim->SetGravity(g);

    // spring constant:
    double C = ui->sfBox->value( ); // spring force constant
    if (sim->springForceConstant != C)
      sim->SetSpringForceConstant(C);

    uint numSteps = ui->frameStepBox->value( );

    time_t t = clock( );

    for (uint i = 0; i < numSteps; i++) {
      if (isRunning == false) {
        numSteps = i;
        break;
      } // if
      this->sim->NextStep(useLLT);
    } // for

    // compute stats:
    double elapsed = double(clock( ) - t) / CLOCKS_PER_SEC;
    double curFPS = double(numSteps) / elapsed;

    if (this->simFPS == 0)
      this->simFPS = curFPS;
    else
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



void MainWindow::on_expButton_clicked( ) {
  this->RunExperiment( );
} // on_expButton_clicked


uint expNum = 0;
double startE;
void MainWindow::RunExperiment( ) {
  double TIJD = 5;
  double minH = 0.01;
  uint numPlots = 5; // how many lines

  if (expNum == 0 && isRunning == false && sim->timeInSim < TIJD) { // first time
    qDebug( ) << "EERSTE KEER !!!";
    this->on_initButton_clicked( );
    sim->SqueezeX(0.95);
    sim->SqueezeY(0.95);

    // gravity:
    bool enableGrav = ui->gravBox->isChecked( );
    double g = (enableGrav)? 0.1 : 0;
    if (sim->gravAcc != g)
      sim->SetGravity(g);

    // spring constant:
    double C = ui->sfBox->value( ); // spring force constant
    if (sim->springForceConstant != C)
      sim->SetSpringForceConstant(C);

    ui->hValSpinBox->setValue(1);
    this->sim->h = ui->hValSpinBox->value( );

    startE = sim->GetTotEnergy( );

    // add new line:
    this->yVecLabels.push_back("#V = " + QString::number(sim->m));
    QVector <double> newYVec;
    this->expVecYs.push_back(newYVec);

    ui->mBox->setValue(ui->mBox->value( ) + 5);

    isRunning = true;
    QTimer::singleShot(10, this, SLOT(RunExperiment( )));
    return;
  } // if

  if (sim->timeInSim >= TIJD) {
    if (isRunning) {
      isRunning = false;
      QTimer::singleShot(250, this, SLOT(RunExperiment( )));
      return;
    } // if

    double curE = sim->GetTotEnergy( );

    this->expVecYs[expVecYs.size( ) - 1].push_back((curE / startE));
    if (expVecYs.size( ) == 1)
      this->expVecX.push_back(log10(ui->hValSpinBox->value( )));

    qDebug( ) << "New datapoint added! :" << expVecX.back( ) << " " << expVecYs[expVecYs.size( ) - 1].back( )
              << ", startE =" << startE << ", curE =" << curE;

    // PLOT:
    QCustomPlot * customPlot = ui->expPlot;

    if (expNum == 0 && expVecYs.size( ) == 1) { // first run, so add title
      customPlot->plotLayout()->insertRow(0);
      customPlot->plotLayout()->addElement(0, 0,
                                           new QCPTextElement(customPlot,
                                           "Energy left after " + QString::number(TIJD) + " seconds, \n#iterations = " + QString::number(sim->iterationsPerStep),
                                                              QFont("sans", 12, QFont::Bold)));
    } // if

    customPlot->clearGraphs( );

    customPlot->setAntialiasedElements(QCP::aePlottables);

    // enable legend and position at top right
    customPlot->legend->setVisible(true);
    customPlot->axisRect( )->insetLayout( )->setInsetAlignment(0, Qt::AlignRight|Qt::AlignTop);


    uint nPlots = this->expVecYs.size( );
    for (uint i = 0; i < nPlots; i++) {
      customPlot->addGraph( );

      customPlot->graph(i)->setData(expVecX, expVecYs[i]);
      customPlot->graph(i)->setPen(this->yVecPens[i]);
      customPlot->graph(i)->setName(this->yVecLabels[i]);
      customPlot->graph(i)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssDisc, 5));
    } // for

    customPlot->xAxis->setLabel("log10 of h [s]");
    customPlot->yAxis->setLabel("Fraction of energy left");
    customPlot->xAxis->setRange(log10(minH), 0);
    customPlot->yAxis->setRange(0, 1);

    customPlot->replot( );

    // reset mesh
    this->on_initButton_clicked( );
    sim->SqueezeX(0.95);
    sim->SqueezeY(0.95);

    // gravity:
    bool enableGrav = ui->gravBox->isChecked( );
    double g = (enableGrav)? 0.1 : 0;
    if (sim->gravAcc != g)
      sim->SetGravity(g);

    // spring constant:
    double C = ui->sfBox->value( ); // spring force constant
    if (sim->springForceConstant != C)
      sim->SetSpringForceConstant(C);

    double lastH = ui->hValSpinBox->value( );
    double newH = ((expNum+1) % 3 == 0)? lastH / 2.5 : lastH / 2;
    ui->hValSpinBox->setValue(newH);

    qDebug( ) << "Changed h value to" << newH << "=" << ui->hValSpinBox->value( );

    startE = sim->GetTotEnergy( );

    if (ui->hValSpinBox->value( ) < minH) { // reset and go to next plot
      qDebug( ) << "< minH, so continue, nPlots =" << nPlots;
      expNum = 0;
      if (nPlots < numPlots)
        QTimer::singleShot(100, this, SLOT(RunExperiment( )));
      return;
    } // if

    isRunning = true;
    expNum++;
    qDebug( ) << "Expnum:" << expNum << ", startE =" << startE;
  } // if

  QWinTaskbarButton *button = new QWinTaskbarButton(this);
  button->setWindow(this->windowHandle( ));
  QWinTaskbarProgress *progress = button->progress( );
  progress->setVisible(true);
  progress->setValue(100 * this->sim->timeInSim / TIJD);

  if (expNum < 4) // act fast, because sim is fast
    QTimer::singleShot(1, this, SLOT(RunExperiment( )));
  else
    QTimer::singleShot(100, this, SLOT(RunExperiment( )));
} // RunExperiment

void MainWindow::UpdateInfoText( ) {
  ui->fpsLabel->setText(QString::number(simFPS, 'f', 0));
  ui->tLabel->setText(QString::number(sim->timeInSim, 'e', 2) + "s");
  ui->netFpsLabel->setText(QString::number(netFPS, 'f', 0));
  ui->simInfoLabel->setText(sim->GetInfoString( ));

  this->UpdateStepPlot( );

  // TEMP!
  QCustomPlot * customPlot = ui->ePlot;

  this->tVec.push_back(sim->timeInSim);
  this->E_grav.push_back(sim->GetGravPotEnergy( ));
  this->E_pot.push_back(sim->GetSpringPotentialEnergy( ));
  this->E_kin.push_back(sim->GetKineticEnergy( ));
  this->E_tot.push_back(E_grav.back( ) + E_kin.back( ) + E_pot.back( ));

  // create graph and assign data to it:
  customPlot->clearGraphs( );

  customPlot->legend->setVisible(true);
  customPlot->axisRect()->insetLayout()->setInsetAlignment(0, Qt::AlignLeft|Qt::AlignTop);

  customPlot->setAntialiasedElements(QCP::aePlottables);

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

  customPlot->graph(0)->setName("Tot E");
  customPlot->graph(1)->setName("E_spring");
  customPlot->graph(2)->setName("E_kin");
  customPlot->graph(3)->setName("E_grav");

  // give the axes some labels:
  customPlot->xAxis->setLabel("Time [s]");
  customPlot->yAxis->setLabel("Energy [J]");

  // set axes ranges, so we see all data:
  customPlot->xAxis->setRange(0, tVec.back( ));
  //customPlot->yAxis->setRange(-fabs(E_tot.back( ) * 1.5), fabs(E_tot.back( ) * 1.5));

  customPlot->graph(0)->rescaleAxes( );
  customPlot->graph(1)->rescaleAxes(true);
  customPlot->graph(2)->rescaleAxes(true);
  customPlot->graph(3)->rescaleAxes(true);

  customPlot->replot( );
} // UpdateInfoText

void MainWindow::UpdateMeshImage( ) {

  bool useAA = ui->aaBox->isChecked( );
  bool drawEdges = ui->drawEdgeBox->isChecked( );
  bool drawSelectedVertices = (ui->MainTab->currentIndex( ) == 2);
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

// updates stepPlot in ui, containing energy loss
// and displacement per step
void MainWindow::UpdateStepPlot( ) {
  if (sim->STEP_disp.size( ) == 0) // nothing to plot
    return;

  if (sim->STEP_disp.size( ) != sim->STEP_totDiffE.size( ))
    qDebug( ) << "Error! STEP vectors not same size!" << sim->STEP_disp.size( ) << sim->STEP_totDiffE.size( );

  QVector <double> stepNrVec;
  for (uint i = 0; i < sim->STEP_disp.size( ); i++)
    stepNrVec.push_back(i + 2);

  //
  // Displacement:
  //
  QCustomPlot * customPlot = ui->stepPlot;
  customPlot->clearGraphs( );
  customPlot->setAntialiasedElements(QCP::aePlottables);

  customPlot->addGraph( );
  customPlot->graph(0)->setData(stepNrVec, QVector<double>::fromStdVector(sim->STEP_disp));
  customPlot->graph(0)->setPen(QPen(Qt::blue));
  customPlot->graph(0)->setName("Displacement");

  customPlot->yAxis->setScaleType(QCPAxis::stLogarithmic);

  customPlot->xAxis->setLabel("Iteration [-]");
  customPlot->yAxis->setLabel("Average displacement [m]");
  customPlot->xAxis->setRange(1, sim->STEP_disp.size( ));
  customPlot->yAxis->setRange(0, 1.5 * sim->STEP_disp[0]);

  //customPlot->graph(0)->rescaleAxes( );
  customPlot->rescaleAxes( );

  customPlot->replot( );

  //
  // Energy:
  //
  customPlot = ui->stepEPlot;
  customPlot->clearGraphs( );
  customPlot->setAntialiasedElements(QCP::aePlottables);

  customPlot->addGraph( );
  customPlot->graph(0)->setData(stepNrVec, QVector<double>::fromStdVector(sim->STEP_totDiffE));
  customPlot->graph(0)->setPen(QPen(Qt::red));
  customPlot->graph(0)->setName("Energy loss");

  customPlot->yAxis->setScaleType(QCPAxis::stLogarithmic);

  customPlot->xAxis->setLabel("Iteration [-]");
  customPlot->yAxis->setLabel("Energy loss [J]");
  customPlot->xAxis->setRange(1, sim->STEP_disp.size( ));
  customPlot->yAxis->setRange(0, 1.5 * sim->STEP_totDiffE[0]);

  customPlot->rescaleAxes( );
  customPlot->replot( );
} // UpdateStepPlot




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

  // gravity:
  bool enableGrav = ui->gravBox->isChecked( );
  double g = (enableGrav)? 0.1 : 0;
  sim->SetGravity(g);

  // spring constant:
  double C = ui->sfBox->value( ); // spring force constant
  sim->SetSpringForceConstant(C);

  this->simFPS = 0;

  on_resetImgButton_clicked( );

  this->tVec.clear( );
  this->E_grav.clear( );
  this->E_pot.clear( );
  this->E_kin.clear( );
  this->E_tot.clear( );

//  this->expVecX.clear( );
//  this->expVecY.clear( );

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


void MainWindow::on_squeezeXButton_clicked( ) {
  this->sim->SqueezeX( );
  UpdateMeshImage( );
} // on_squeezeXButton_clicked


void MainWindow::on_squeezeYButton_clicked( ) {
  this->sim->SqueezeY( );
  UpdateMeshImage( );
} // on_squeezeYButton_clicked

