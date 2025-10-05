// Fill out your copyright notice in the Description page of Project Settings.


#include "POIDataAsset.h"

void UPOIDataAsset::PostEditChangeProperty(FPropertyChangedEvent& PropertyChangedEvent)
{
	Super::PostEditChangeProperty(PropertyChangedEvent);

	for (FPOI& POI : POIList) {
		if (POI.ID.IsNone()) {
			POI.ID = FName(*FGuid::NewGuid().ToString());
		}
	}
}