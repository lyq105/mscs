#ifndef ANALYSISTYPE_H
#define ANALYSISTYPE_H

#include "ui_analysistype.h"

class QString;
class AnalysisType : public QDialog
{
    Q_OBJECT

public:
    explicit AnalysisType(QWidget *parent = 0);
    QString get_prj_name();
    QString get_analysis_type();
    QString get_prj_folder();

protected:
    void changeEvent(QEvent *e);

private slots:
    void open_prj_folder();

private:
    Ui::AnalysisType ui;
};

#endif // ANALYSISTYPE_H
