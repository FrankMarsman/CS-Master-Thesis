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

#include <QClipboard>
#include <QPixmap>
#include <QFileDialog>

using namespace Eigen;
using namespace std;

MainWindow::MainWindow(QWidget *parent) :
  QMainWindow(parent),
  ui(new Ui::MainWindow) {
  ui->setupUi(this);

  uint N = 30;
  double l = ui->meshSizeBox->value( );
  double mass = ui->meshMassBox->value( );
  this->sim = new Sim2D(N * N, mass, l);

  ui->mBox->setValue(N);

  // make sure ui updates when settings change:
  QObject::connect(ui->aaBox, SIGNAL(toggled(bool)), SLOT(UpdateMeshImage( )));
  QObject::connect(ui->drawEdgeBox, SIGNAL(toggled(bool)), SLOT(UpdateMeshImage( )));
  QObject::connect(ui->drawModeBox, SIGNAL(currentIndexChanged(int)), SLOT(UpdateMeshImage( )));

  QObject::connect(ui->sX1Slider, SIGNAL(actionTriggered(int)), SLOT(UpdateSelectedVertices( )));
  QObject::connect(ui->sX2Slider, SIGNAL(actionTriggered(int)), SLOT(UpdateSelectedVertices( )));
  QObject::connect(ui->sY1Slider, SIGNAL(actionTriggered(int)), SLOT(UpdateSelectedVertices( )));
  QObject::connect(ui->sY2Slider, SIGNAL(actionTriggered(int)), SLOT(UpdateSelectedVertices( )));
  // floor related controls:
  QObject::connect(ui->ffDistBox, SIGNAL(valueChanged(double)), SLOT(UpdateFloor( )));
  QObject::connect(ui->ffPosBox, SIGNAL(valueChanged(double)), SLOT(UpdateFloor( )));
  QObject::connect(ui->ffStrengthBox, SIGNAL(valueChanged(double)), SLOT(UpdateFloor( )));
  QObject::connect(ui->ffCheckBox, SIGNAL(toggled(bool)), SLOT(UpdateFloor( )));
  // view related controls:
  QObject::connect(ui->imgXBox, SIGNAL(valueChanged(double)), SLOT(UpdateImgParams( )));
  QObject::connect(ui->imgYBox, SIGNAL(valueChanged(double)), SLOT(UpdateImgParams( )));
  QObject::connect(ui->imgWBox, SIGNAL(valueChanged(double)), SLOT(UpdateImgParams( )));

  QObject::connect(ui->imgScaleBox, SIGNAL(valueChanged(double)), SLOT(UpdateMeshImageInfo( )));



  prevRunStep = 0;
  simStepTimer = clock( );

  SIM_MAX_RUNTIME = 1e99; // default no max time

  // init pens for plots:
  this->yVecPens.push_back(QPen(Qt::blue));
  this->yVecPens.push_back(QPen(Qt::red));
  this->yVecPens.push_back(QPen(Qt::green));
  this->yVecPens.push_back(QPen(Qt::yellow));
  this->yVecPens.push_back(QPen(Qt::magenta));


  isRunning = true; // set to opposite what you want

  on_nextButton_clicked( );
  on_resetImgButton_clicked( );
  on_meshTypeBox_currentIndexChanged(ui->meshTypeBox->currentIndex( ));

  UpdateSelectedVertices( );
  UpdateMeshImage( );

  QTimer::singleShot(10, this, SLOT(on_initButton_clicked( )));
  QTimer::singleShot(250, this, SLOT(Run( )));

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
  double g = (enableGrav)? ui->gValBox->value( ) : 0;
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
  this->sim->enableFriction = ui->frickCheckBox->isChecked( );
  this->sim->fricC = ui->fricSpinBox->value( );

  bool useLLT = ui->lltBox->isChecked( );

  bool useTimeLimit = ui->timeLimitCheckBox->isChecked( );
  if (useTimeLimit)
    SIM_MAX_RUNTIME = ui->maxTimeBox->value( );

  if (isRunning) {

    // gravity:
    bool enableGrav = ui->gravBox->isChecked( );
    double g = (enableGrav)? ui->gValBox->value( ) : 0;
    if (sim->gravAcc != g)
      sim->SetGravity(g);

    // spring constant:
    double C = ui->sfBox->value( ); // spring force constant
    if (sim->springForceConstant != C)
      sim->SetSpringForceConstant(C);

    uint numSteps = ui->frameStepBox->value( );

    time_t t = clock( );

    for (uint i = 0; i < numSteps; i++) {
      if (isRunning == false || (useTimeLimit && sim->timeInSim >= SIM_MAX_RUNTIME)) {
        numSteps = i;
        // set isrunning to false:
        isRunning = true;
        on_nextButton_clicked( ); // this makes isRunning false
        break;
      } // if
      else
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
    if (isnan(simFPS))
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
  SIM_MAX_RUNTIME = 1;
  ui->maxTimeBox->setValue(SIM_MAX_RUNTIME);
  ui->timeLimitCheckBox->setChecked(true);
  this->RunExperiment( );
} // on_expButton_clicked


uint expNum = 0;
double startE;
void MainWindow::RunExperiment( ) {
  double minH = 0.0001;
  uint numPlots = 5; // how many lines

  QWinTaskbarButton *button = new QWinTaskbarButton(this);
  button->setWindow(this->windowHandle( ));
  QWinTaskbarProgress *progress = button->progress( );

  if (expNum == 0 && isRunning == false && sim->timeInSim < SIM_MAX_RUNTIME) { // first time
    if (expVecYs.size( ) > 0) { // not REAL first time
      //ui->mBox->setValue(ui->mBox->value( ) * 2);
      ui->meshMassBox->setValue(ui->meshMassBox->value( ) * 4);
    } // if
    this->on_initButton_clicked( );

    qDebug( ) << "EERSTE KEER !!! m =" << sim->m;

    ui->hValSpinBox->setValue(0.1);
    this->sim->h = ui->hValSpinBox->value( );

    startE = sim->GetTotEnergy( );

    // add new line:
    //this->yVecLabels.push_back("#V = " + QString::number(sim->m) + ", E₀ = " + QString::number(startE, 'g', 2) + "J");
    this->yVecLabels.push_back("M = " + QString::number(ui->meshMassBox->value( )) + "kg, E₀ = " + QString::number(startE, 'g', 2) + "J");
    QVector <double> newYVec;

    if (expVecYs.size( ) > 0) { // not first vector, so make same size
      for (int i = 0; i < expVecYs[0].size( ); i++)
        newYVec.push_back(0);
    } // if

    this->expVecYs.push_back(newYVec);

    progress->setVisible(true);
    progress->setValue(0);

    isRunning = true;
    QTimer::singleShot(10, this, SLOT(RunExperiment( )));
    return;
  } // if

  if (sim->timeInSim >= SIM_MAX_RUNTIME) {
    if (isRunning) {
      isRunning = false;
      QTimer::singleShot(250, this, SLOT(RunExperiment( )));
      return;
    } // if

    double curE = sim->GetTotEnergy( );

    if (expVecYs.size( ) == 1) {
      this->expVecX.push_back(-log10(ui->hValSpinBox->value( )));
      this->expVecYs[expVecYs.size( ) - 1].push_back((curE / startE));
    } // if
    else { // yvec is already filled, so set value instead of add
      this->expVecYs[expVecYs.size( ) - 1][expNum] = ((curE / startE));
    } // else

    qDebug( ) << "New datapoint added! :" << expVecX.back( ) << " " << expVecYs[expVecYs.size( ) - 1].back( )
              << ", startE =" << startE << ", curE =" << curE << "simTime =" << sim->timeInSim;

    // PLOT:
    QCustomPlot * customPlot = ui->expPlot;

    if (expNum == 0 && expVecYs.size( ) == 1) { // first run, so add title
      customPlot->plotLayout()->insertRow(0);
      customPlot->plotLayout()->addElement(0, 0,
                                           new QCPTextElement(customPlot,
                                           "Energy left after " + QString::number(SIM_MAX_RUNTIME) + " seconds, \n#iterations = " + QString::number(sim->iterationsPerStep) + ", #V = " + QString::number(sim->m),
                                                              QFont("sans", 12, QFont::Bold)));
    } // if

    customPlot->clearGraphs( );

    customPlot->setAntialiasedElements(QCP::aePlottables);

    // enable legend and position at bottom right
    customPlot->legend->setVisible(true);
    customPlot->axisRect( )->insetLayout( )->setInsetAlignment(0, Qt::AlignLeft|Qt::AlignTop);


    uint nPlots = this->expVecYs.size( );
    for (uint i = 0; i < nPlots; i++) {
      customPlot->addGraph( );

      customPlot->graph(i)->setData(expVecX, expVecYs[i]);
      customPlot->graph(i)->setPen(this->yVecPens[i]);
      customPlot->graph(i)->setName(this->yVecLabels[i]);
      customPlot->graph(i)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssDisc, 5));
    } // for

    customPlot->xAxis->setLabel("log10(1/h) [s]");
    customPlot->yAxis->setLabel("Fraction of energy left");
    customPlot->xAxis->setRange(0, -log10(minH));
    customPlot->yAxis->setRange(0, 1);

    customPlot->replot( );

    // reset mesh
    this->on_initButton_clicked( );

    double lastH = ui->hValSpinBox->value( );
    double newH = ((expNum+1) % 3 == 0)? lastH / 2.5 : lastH / 2;
    //double newH = (expNum % 2 == 0)? lastH / 2 : lastH / 5;
    ui->hValSpinBox->setValue(newH);

    qDebug( ) << "Changed h value to" << newH << "=" << ui->hValSpinBox->value( );

    startE = sim->GetTotEnergy( );

    if (ui->hValSpinBox->value( ) < minH) { // reset and go to next plot
      qDebug( ) << "< minH, so continue, nPlots =" << nPlots;
      expNum = 0;
      progress->setVisible(false);
      progress->setValue(0);
      if (nPlots < numPlots)
        QTimer::singleShot(100, this, SLOT(RunExperiment( )));
      return;
    } // if

    isRunning = true;
    expNum++;
    qDebug( ) << "Expnum:" << expNum << ", startE =" << startE;
  } // if

  progress->setVisible(true);
  progress->setValue(100 * this->sim->timeInSim / SIM_MAX_RUNTIME);

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

  bool plotThings = ui->plotCheckBox->isChecked( );

  if (plotThings) {
  this->UpdateStepPlot( );
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

    customPlot->graph(0)->rescaleAxes( );
    customPlot->graph(1)->rescaleAxes(true);
    customPlot->graph(2)->rescaleAxes(true);
    customPlot->graph(3)->rescaleAxes(true);

    customPlot->replot( );
  } // if
} // UpdateInfoText

void MainWindow::UpdateMeshImage( ) {

  bool useAA = ui->aaBox->isChecked( );
  bool drawEdges = ui->drawEdgeBox->isChecked( );
  bool drawSelectedVertices = (ui->MainTab->currentIndex( ) == 2);
  bool drawLockedVertices = (ui->MainTab->currentIndex( ) == 0);

  uint DRAW_MODE = ui->drawModeBox->currentIndex( );


  //qDebug( ) << "UpdateMeshImage, DRAW_MODE =" << DRAW_MODE;

  QImage img = sim->ToQImage(700, useAA, drawEdges, drawSelectedVertices, drawLockedVertices, DRAW_MODE);
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

// updates floor parameters and redraws mesh
void MainWindow::UpdateFloor(bool redraw) {
  bool enableFloor = ui->ffCheckBox->isChecked( );

  double fPos = ui->ffPosBox->value( );
  double fDist = ui->ffDistBox->value( );
  double fStrength = ui->ffStrengthBox->value( );

  sim->floorEnabled = enableFloor;
  sim->floorC = fStrength;
  sim->floorDist = fDist;
  sim->floorHeight = fPos;
  sim->gravRefHeight = fPos;

  if (redraw)
    UpdateMeshImage( );
} // UpdateFloor


// updates center of mesh image and size of view
void MainWindow::UpdateImgParams( ) {
  double viewSize = ui->imgWBox->value( );
  double cx = ui->imgXBox->value( );
  double cy = ui->imgYBox->value( );

  sim->imgCenterX = cx;
  sim->imgCenterY = cy;
  sim->imgViewSize = viewSize;

  UpdateMeshImage( );
} // UpdateImgParams



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
  double l = ui->meshSizeBox->value( );
  double mass = ui->meshMassBox->value( );
  bool useDSprings = ui->diagSpringCheckBox->isChecked( );
  m *= m; // make squared

  delete this->sim;
  this->sim = NULL;

  uint mType = ui->meshTypeBox->currentIndex( );

  if (mType == 0) { // square mesh
    this->sim = new Sim2D(m, mass, l, useDSprings);
  } // if

  if (mType == 1) { // mesh from image
    QImage img;
    QString str = ui->meshFileEdit->text( );
    if (str.length( ) == 0) {
      qDebug( ) << "Error! No file selected";
      return;
    } // if
    bool r = img.load(str);
    if (r == false) {
      qDebug( ) << "Error! Image load failed";
      return;
    } // if
    double sFac = ui->imgScaleBox->value( );
    img = img.scaled(img.width( ) * sFac, img.height( ) * sFac);
    this->sim = new Sim2D(img, mass, l, useDSprings);
    ui->statusBar->showMessage("Mesh from image, m = " + QString::number(sim->m));
  } // if

  if (mType == 2) { // spring
    this->sim = new Sim2D_Spring(10, 50, 0.1, 1, mass);
  } // if

  PerformInitDeform( );

  // gravity:
  bool enableGrav = ui->gravBox->isChecked( );
  double g = (enableGrav)? ui->gValBox->value( ) : 0;
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

  UpdateFloor(false);
  UpdateSelectedVertices( );
  UpdateMeshImage( );

  ui->statusBar->showMessage("Done initializing mesh. m = " + QString::number(sim->m)
                             + ", #springs = " + QString::number(sim->numEdges)
                             + ", E = " + QString::number(sim->GetTotEnergy( )) + "J");

  ui->simInfoLabel->setText(sim->GetInfoString( ));
} // on_initButton_clicked

// checks whether initial set of deformations
// has to be done (ui.initActionBox)
void MainWindow::PerformInitDeform( ) {
  uint actionIndex = ui->initActionBox->currentIndex( );
  uint N = sqrt(sim->m);
  uint m = sim->m;
  uint teller;

  if (actionIndex == 0) // do nothing
    return;

  if (actionIndex == 1) { // lock bottom vertices
    sim->lockedVertices.clear( );
    for (uint i = 0; i < N; i++) // lock first N vertices
      sim->lockedVertices.push_back(i);
    sim->gravRefHeight += 1;
  } // if

  if (actionIndex == 2) { // add opposite x-velocity
    teller = N * int(N / 4);
    // add velocity to right for top:
    sim->selectedVertices.clear( );
    for (uint i = 0; i < teller; i++) // select top 25%
      sim->selectedVertices.push_back(m - i - 1);
    sim->AddVelocity(0.2, 0);
    // add velocity to left for bottom:
    sim->selectedVertices.clear( );
    for (uint i = 0; i < teller; i++) // select top 25%
      sim->selectedVertices.push_back(i);
    sim->AddVelocity(-0.2, 0);
  } // if

  if (actionIndex == 3) { // add opposite x-velocity and opposite y-velocity
    double speed = 0.2;
    sim->SetSelectedVertices(0, 1, 0, 0.25);
    sim->AddVelocity(speed, 0);
    sim->SetSelectedVertices(0, 1, 0.75, 1);
    sim->AddVelocity(-speed, 0);
    sim->SetSelectedVertices(0, 0.25, 0, 1);
    sim->AddVelocity(0, speed);
    sim->SetSelectedVertices(0.75, 1, 0, 1);
    sim->AddVelocity(0, -speed);
  } // if

  if (actionIndex == 4) { // same as 3, but with x and y squeeze
    double speed = 0.2;
    sim->SetSelectedVertices(0, 1, 0, 0.25);
    sim->AddVelocity(speed, 0);
    sim->SetSelectedVertices(0, 1, 0.75, 1);
    sim->AddVelocity(-speed, 0);
    sim->SetSelectedVertices(0, 0.25, 0, 1);
    sim->AddVelocity(0, speed);
    sim->SetSelectedVertices(0.75, 1, 0, 1);
    sim->AddVelocity(0, -speed);
    sim->SqueezeX(0.9);
    sim->SqueezeY(0.9);
  } // if

  if (actionIndex == 5) { // lock bottom and top vertices and strech out
    sim->lockedVertices.clear( );
    for (uint i = 0; i < N; i++) { // lock first and last N vertices
      sim->lockedVertices.push_back(i);
      sim->lockedVertices.push_back(m - i - 1);
    } // for
    sim->SqueezeY(1.5);
  } // if
} // PerformInitDeform

// user selects different tab
void MainWindow::on_MainTab_currentChanged(int index) {
  if (index < 0)
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
  ui->simInfoLabel->setText(sim->GetInfoString( ));
} // on_squeezeXButton_clicked


void MainWindow::on_squeezeYButton_clicked( ) {
  this->sim->SqueezeY( );
  UpdateMeshImage( );
  ui->simInfoLabel->setText(sim->GetInfoString( ));
} // on_squeezeYButton_clicked

// copies content of expPlot to clipboard
void MainWindow::on_copyExpBut_clicked( ) {
  QPixmap pm = ui->expPlot->toPixmap(600, 450);
  QApplication::clipboard( )->setPixmap(pm);
} // on_copyExpBut_clicked


// copies mesh image to clipboard
void MainWindow::on_copyMeshImgButton_clicked( ) {
  const QPixmap* pm = ui->imgLabel->pixmap( );
  QApplication::clipboard( )->setPixmap(*pm);
} // on_copyMeshImgButton_clicked

// mesh type changed, if mesh type is from image, then show
// image selection box
void MainWindow::on_meshTypeBox_currentIndexChanged(int index) {
  if (index == 1)
    ui->meshFileGroupBox->setVisible(true);
  else
    ui->meshFileGroupBox->setVisible(false);
} // on_meshTypeBox_currentIndexChanged

// opens dialog that let's user choose mesh file
void MainWindow::on_chooseMeshFileButton_clicked( ) {
  QString str = QFileDialog::getOpenFileName(
        this,
        tr("Open mesh image"), "",
        tr("Image (*.*);;All Files (*)"));
  ui->meshFileEdit->setText(str);

  UpdateMeshImageInfo( );
} // on_chooseMeshFileButton_clicked

// computes info about mesh image
void MainWindow::UpdateMeshImageInfo( ) {
  QImage img;
  QString str = ui->meshFileEdit->text( );
  if (str.length( ) == 0) {
    ui->statusBar->showMessage("Error! No file selected");
    return;
  } // if
  bool r = img.load(str);
  if (r == false) {
    ui->statusBar->showMessage("Error! Image load failed");
    return;
  } // if
  double sFac = ui->imgScaleBox->value( );
  img = img.scaled(img.width( ) * sFac, img.height( ) * sFac);
  uint imgW = img.width( ), imgH = img.height( );
  QColor col;
  uint m = 0;
  // determine number of vertices [m]:
  for (uint x = 1; x < imgW - 1; x++)
    for (uint y = 1; y < imgH - 1; y++) {
      col = img.pixel(x, y);
      if (col.red( ) > 0) // non-black pixel
        m++;
    } // for
  ui->statusBar->showMessage("Selected image contains " + QString::number(m) + " vertices");
} // UpdateMeshImageInfo

// shows image of constant system matrix
void MainWindow::on_updateLMatImg_clicked( ) {
  QImage img = sim->GetlMatrixImage( );
  ui->lMatImgLabel->setPixmap(QPixmap::fromImage(img).scaled(500, 500));
} // on_updateLMatImg_clicked





