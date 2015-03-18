//
//    PEDPArguments.c
//
//    soap2-dp
//
//    Copyright (C) 2013, HKU
//
//    This program is free software; you can redistribute it and/or
//    modify it under the terms of the GNU General Public License
//    as published by the Free Software Foundation; either version 2
//    of the License, or (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program; if not, write to the Free Software
//    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
//
///////////////////////////////////////////////////////////////////////////////////////////////
#include <stdlib.h>

#include "PEDPArguments.h"

PEDPSetting * PEDPSettingConstruct () {
    PEDPSetting * pedpSettings = (PEDPSetting*) malloc (sizeof(PEDPSetting));
    
    pedpSettings->SGAScoreTF = 0.30f;
    pedpSettings->SGASoftHeadClipLength = 0.0f;
    pedpSettings->SGASoftTailClipLength = 0.40f;
    pedpSettings->SGASoftTotalClipLength = 0.40f;
    pedpSettings->SGALongReadSingleEndDPScoreTF = 0.30f;
    pedpSettings->SGALongReadDefaultDPScoreTF = 0.30f;
    
    pedpSettings->SGAOrphanEnhancement = 0;
    pedpSettings->SGAOrphanTriggerTF = 0.30f;

    pedpSettings->SGAOrphanExtendEnhancement = 0;
    pedpSettings->SGAOrphanExtendTriggerTF = 0.30f;
    
    pedpSettings->SGASeedEnhancement_Sweep = 0;
    pedpSettings->SGASeedLength_Sweep = 80.0f;
    pedpSettings->SGASeedLengthOverlap_Sweep = 0.15f;
    pedpSettings->SGASeedLooseCriteria_Sweep = -1;
    pedpSettings->SGASeedOccLimitation_Sweep = 4096;
    
    pedpSettings->SGASeedEnhancement = 0;
    pedpSettings->SGASeedLength = 0.30f;
    pedpSettings->SGASeedLengthOverlap = 0.15f;
    pedpSettings->SGASeedLooseCriteria = -1;
    pedpSettings->SGASeedOccLimitation = 4096;
    pedpSettings->SGASeedLengthDPThreshold_Sweep = 1;
    
    pedpSettings->SGASeedEnhancement_1 = 0;
    pedpSettings->SGASeedLength_1 = 0.30f;
    pedpSettings->SGASeedLengthOverlap_1 = 0.15f;
    pedpSettings->SGASeedLooseCriteria_1 = -1;
    pedpSettings->SGASeedOccLimitation_1 = 4096;
    
    pedpSettings->SGASeedEnhancement_2 = 0;
    pedpSettings->SGASeedLength_2 = 0.30f;
    pedpSettings->SGASeedLengthOverlap_2 = 0.15f;
    pedpSettings->SGASeedLooseCriteria_2 = -1;
    pedpSettings->SGASeedOccLimitation_2 = 4096;
    return pedpSettings;
}

void PEDPSettingFree(PEDPSetting * pedpSettings) {
    free(pedpSettings);
}
