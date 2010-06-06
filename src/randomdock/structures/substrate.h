/**********************************************************************
  RandomDock - Substrate: Wrapper for Avogadro::Molecule to hold the 
  central molecule in a docking problem

  Copyright (C) 2009 by David C. Lonie

  This file is part of the Avogadro molecular editor project.
  For more information, see <http://avogadro.openmolecules.net/>

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation version 2 of the License.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
 ***********************************************************************/

#ifndef SUBSTRATEMOL_H
#define SUBSTRATEMOL_H

#include <globalsearch/structure.h>

#include <QDebug>
#include <QDateTime>
#include <QTextStream>

using namespace GlobalSearch;

namespace RandomDock {

  class Substrate : public Structure
  {
    Q_OBJECT

   public:
    Substrate(QObject *parent = 0);
    Substrate(Molecule *mol);
    virtual ~Substrate();

   signals:

   public slots:
    double prob(uint index) { checkProbs(); return m_probs.at(index);};
    void sortConformers();
    void generateProbabilities();
    void checkProbs() {if ((uint)m_probs.size() != numConformers()) generateProbabilities();};
    int getRandomConformerIndex();

   private slots:

   private:
    QList<double> m_probs;
  };

} // end namespace Avogadro

#endif