#include "analysistype.h"
#include <QFileDialog>
AnalysisType::AnalysisType(QWidget *parent) :
    QDialog(parent)
{
    ui.setupUi(this);
    connect(ui.tb_open_prj_folder,SIGNAL(clicked()), this, SLOT(open_prj_folder()));
}

void AnalysisType::changeEvent(QEvent *e)
{
    QDialog::changeEvent(e);
    switch (e->type()) {
    case QEvent::LanguageChange:
        ui.retranslateUi(this);
        break;
    default:
        break;
    }
}


QString AnalysisType::get_prj_name()
{
    return ui.lineEdit_prj_name->text();
}
QString AnalysisType::get_analysis_type()
{
    return "elastic";
}
QString AnalysisType::get_prj_folder()
{
    return ui.lineEdit_folder->text();
}

void AnalysisType::open_prj_folder()
{
    QFileDialog* openFilePath = new QFileDialog( this, tr("Please Choose Folder"), "file");
    openFilePath->setFileMode( QFileDialog::DirectoryOnly);
    if ( openFilePath->exec() == QDialog::Accepted )
    {
        ui.lineEdit_folder->setText( openFilePath->selectedFiles()[0]);
    }
    delete openFilePath;
}
