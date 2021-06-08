#-------------------------------------------------
#
# Project created by QtCreator 2021-01-27T19:44:40
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets printsupport

TARGET = SFRTest
TEMPLATE = app

# The following define makes your compiler emit warnings if you use
# any feature of Qt which has been marked as deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.
DEFINES += QT_DEPRECATED_WARNINGS

# You can also make your code fail to compile if you use deprecated APIs.
# In order to do so, uncomment the following line.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

CONFIG += c++11

SOURCES += \
        main.cpp \
        mainwindow.cpp \
    sfr.cpp \
    qcustomplot.cpp

HEADERS += \
        mainwindow.h \
    sfr.h \
    qcustomplot.h

FORMS += \
        mainwindow.ui
#LIBS += -L"$$PWD/lib" -lopencv_world440
#LIBS += -L"$$PWD/lib" -lopencv_world440d
INCLUDEPATH += $$PWD/include

CONFIG(debug, debug|release): {
LIBS += -L$$PWD/lib \
-lopencv_core249d \
-lopencv_imgproc249d \
-lopencv_highgui249d \
-lopencv_ml249d \
-lopencv_video249d \
-lopencv_features2d249d \
-lopencv_calib3d249d \
-lopencv_objdetect249d \
-lopencv_contrib249d \
-lopencv_legacy249d \
-lopencv_flann249d \
-lopencv_photo249d \
-lopencv_ocl249d \
-lopencv_imgproc249d\
-lopencv_legacy249d
} else:CONFIG(release, debug|release): {
LIBS += -L$$PWD/lib \
-lopencv_core249 \
-lopencv_imgproc249 \
-lopencv_highgui249 \
-lopencv_ml249 \
-lopencv_video249 \
-lopencv_features2d249 \
-lopencv_calib3d249 \
-lopencv_objdetect249 \
-lopencv_contrib249 \
-lopencv_legacy249 \
-lopencv_flann249 \
-lopencv_photo249\
-lopencv_ocl249\
-lopencv_imgproc249 \
-lopencv_legacy249
}
# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target
