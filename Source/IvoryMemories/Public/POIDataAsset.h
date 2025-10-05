// Fill out your copyright notice in the Description page of Project Settings.

#pragma once

#include "CoreMinimal.h"
#include "Engine/DataAsset.h"
#include "POIStruct.h"
#include "POIDataAsset.generated.h"

UCLASS(BlueprintType)
class IVORYMEMORIES_API UPOIDataAsset : public UDataAsset
{
	GENERATED_BODY()

public:
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "POI List")
	TArray<FPOI> POIList;

	virtual void PostEditChangeProperty(FPropertyChangedEvent& PropertyChangedEvent) override;
	
};
