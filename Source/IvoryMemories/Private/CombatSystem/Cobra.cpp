// Fill out your copyright notice in the Description page of Project Settings.


#include "CombatSystem/Cobra.h"
#include "CombatSystem/Cell.h"

ACobra::ACobra()
{
	MaxMoveDistance = 3;
	Rank = 3;
	GenerateUniqueName(TEXT("Cobra"), 1);
	// HERE set the mesh to the unit
}


bool ACobra::CanMoveTo(ACell* TargetCell) const
{
    if (!Super::CanMoveTo(TargetCell)) { return false; }

    // specific logic: diagonal move
    if (CurrentCell) {
        FVector2D Delta = TargetCell->Coordinates - CurrentCell->Coordinates;
        return FMath::Abs(Delta.X) == FMath::Abs(Delta.Y); // diagonal only
    }
    return true;
}

void ACobra::PoisonAbility(AFigure* TargetFigure)
{
    UE_LOG(LogTemp, Warning, TEXT("Cobra poisoned target"));
}