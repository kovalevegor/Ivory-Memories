// Fill out your copyright notice in the Description page of Project Settings.

#pragma once

#include "CoreMinimal.h"
#include "GameFramework/Actor.h"
#include "Figure.generated.h"

class ACell;

UCLASS()
class IVORYMEMORIES_API AFigure : public AActor
{
	GENERATED_BODY()
	
public:	
	// Sets default values for this actor's properties
	AFigure();

protected:
	// Called when the game starts or when spawned
	//virtual void BeginPlay() override;

public:	
	// Called every frame
	//virtual void Tick(float DeltaTime) override;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Figure")
	ACell* CurrentCell; // cell reference for bijection

	//------------------------------------------------------------

	UFUNCTION(BlueprintCallable, Category = "Figure")
	void SetCell(ACell* newCell); 

};
