// Fill out your copyright notice in the Description page of Project Settings.

#pragma once

#include "CoreMinimal.h"
#include "CombatSystem/Figure.h"
#include "Cobra.generated.h"


// UNIT Cobra class INHERITS AFigure class 

UCLASS()
class IVORYMEMORIES_API ACobra : public AFigure
{
	GENERATED_BODY()
	
public:
	ACobra();

	// ---INHERITED---
	virtual bool CanMoveTo(ACell* TargerCell) const override;
	
	// ---PROPERTIES---
	UFUNCTION(BlueprintCallable, Category = "Cobra")
	void PoisonAbility(AFigure* TargetFigure);
};
