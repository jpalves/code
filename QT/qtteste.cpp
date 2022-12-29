#include <QCoreApplication>
#include <QApplication>
#include <QPushButton>

int main(int argc, char** argv){
  QApplication app(argc, argv);

  QPushButton hello("Hello world!", 0);
  hello.resize(800, 600);
  hello.show();

  return app.exec();
}
