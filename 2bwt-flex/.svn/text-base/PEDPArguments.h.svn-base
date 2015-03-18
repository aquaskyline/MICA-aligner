//
//    PEDPArguments.h
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

#ifndef __PEDPARGUMENTS_H__
#define __PEDPARGUMENTS_H__

typedef struct PEDPSetting {
    double SGAScoreTF;
    double SGASoftHeadClipLength;
    double SGASoftTailClipLength;
    double SGASoftTotalClipLength;
    double SGALongReadSingleEndDPScoreTF;
    double SGALongReadDefaultDPScoreTF;
    
    int    SGAOrphanEnhancement;
    double SGAOrphanTriggerTF;

    int    SGAOrphanExtendEnhancement;
    double SGAOrphanExtendTriggerTF;

    int    SGASeedEnhancement_Sweep;
    double SGASeedLength_Sweep;
    double SGASeedLengthOverlap_Sweep;
    int    SGASeedLooseCriteria_Sweep;
    int    SGASeedOccLimitation_Sweep;
    double SGASeedLengthDPThreshold_Sweep;

    int    SGASeedEnhancement;
    double SGASeedLength;
    double SGASeedLengthOverlap;
    int    SGASeedLooseCriteria;
    int    SGASeedOccLimitation;

    int    SGASeedEnhancement_1;
    double SGASeedLength_1;
    double SGASeedLengthOverlap_1;
    int    SGASeedLooseCriteria_1;
    int    SGASeedOccLimitation_1;

    int    SGASeedEnhancement_2;
    double SGASeedLength_2;
    double SGASeedLengthOverlap_2;
    int    SGASeedLooseCriteria_2;
    int    SGASeedOccLimitation_2;
} PEDPSetting;

PEDPSetting * PEDPSettingConstruct();
void PEDPSettingFree(PEDPSetting * pedpSettings);

#endif

