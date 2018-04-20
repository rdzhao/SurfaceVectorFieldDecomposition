#pragma once

#include <QWidget>
#include <QMainWindow>
#include <QOpenGLWidget>
#include <QFileDialog>
#include <QMessageBox>

#include <QMenuBar>
#include <QStatusBar>
#include <QToolBar>
#include <QIcon>

#include <QTextStream>

#include <iostream>

#include "MyOpenGLWidget.h"

using namespace std;

class Viewer : public QMainWindow
{
	Q_OBJECT

public:
	Viewer(QWidget* parent = Q_NULLPTR);

	private slots:
	void load();

	void showFaceSlot();
	void showEdgeSlot();

	void faceSelectionSlot();
	void edgeSelectionSlot();
	void vertSelectionSlot();

	void showVectorFieldSlot();
	void showExactSlot() {};
	void showCoexactSlot() {};
	void showHarmonicSlot() {};

	void clearSceneSlot();

	void decompisitionSlot();

private:
	void createCanvas();
	void createMenu();
	void createStatus();
	void createToolBar();

	//widget
private:
	OGLWidget* canvas;

	// menu
private:
	QMenu *fileMenu;
	QMenu *sceneMenu;
	QMenu *selectMenu;
	QMenu *opMenu;

	QAction *newAct;
	QAction *openAct;
	QAction *saveAct;
	QAction *exitAct;

	QAction *showOriginalFieldAct;
	QAction *showExactAct;
	QAction *showCoexactAct;
	QAction *showHarmonicAct;
	QAction *ClearAct;

	QAction *decomAct;

private:
	QToolBar * toolBar;

	QAction *tbImportAct;
	QAction *tbExportAct;

	QAction *showFaceAct;
	QAction *showEdgeAct;

	QAction *tbSelectFaceAct;
	QAction *tbSelectEdgeAct;
	QAction *tbSelectVertAct;
	
	//data
private:
	Mesh* mesh;
	Mesh* ball;
};