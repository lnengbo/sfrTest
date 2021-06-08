#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "qcustomplot.h"
MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    QCustomPlot *custplot = new QCustomPlot();
}

MainWindow::~MainWindow()
{
    delete ui;
}
