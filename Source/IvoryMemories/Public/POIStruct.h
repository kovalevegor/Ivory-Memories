// Fill out your copyright notice in the Description page of Project Settings.

#pragma once

#include "CoreMinimal.h"
#include "POIStruct.generated.h"

USTRUCT(BlueprintType)
struct FPOI
{
	GENERATED_BODY()

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "POI")
	FName ID;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "POI")
	FString Name;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "POI")
	int32 ClusterID;

};

//**
// * 
// */
//class IVORYMEMORIES_API POIStruct
//{
//public:
//	POIStruct();
//	~POIStruct();
//};
