/* Copyright (c) 2009, NICTA
 * All rights reserved.
 *
 * The contents of this file are subject to the Mozilla Public License
 * Version 1.1 (the "License"); you may not use this file except in
 * compliance with the License. You may obtain a copy of the License at
 * http://www.mozilla.org/MPL/
 *
 * Software distributed under the License is distributed on an "AS IS"
 * basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See the
 * License for the specific language governing rights and limitations
 * under the License.
 *
 * Authors: Jin Yu (jin.yu@nicta.com.au)
 *          Choon Hui Teo (choonhui.teo@anu.edu.au)
 *
 * Created: (14/02/2008)
 *
 * Last Updated: 16/01/2009
 */

#ifndef _LINESEARCH_CPP_
#define _LINESEARCH_CPP_

#include "linesearch.hpp"

CLinesearch::CLinesearch(CModel* model, CLoss* loss, CData* data)
{
   _model = model;
   numElements = data->slice_size();
   C = loss->GetScalingFactor();
   numClass = model->GetNumOfW();
   eta = 1.0;
   error = 0.0;
}


#endif
