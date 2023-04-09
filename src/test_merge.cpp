void segMerge(int minSize, int minDiff){

    // Merge de l'image segmenté en fonction de la distance entre leurs moyennes
    vector<region> smallDist;
    // Récupération de la liste des régions dont la distance est inférieur à la distance minimum
    float diff;
    for(region:listRegion){
        if(diff(region, voisin0)<minDiff && region.size!=0)
            if(region.size<voisin0.size) smallDist.push_back(region);
            else smallDist.push_back(voisin0);
    }
    sort(smallDist);
    
    // Merge de la liste des régions
    while(smallDist.size!=0){
        mergeRegion(currentRegion, voisin0);
        update(smallDist, voisin0);
    }



    // Merge de l'image segmenté en fonction de la taille des régions
    vector<region> listSmall;
    // Récupération de la liste des régions inférieures à la taille minimum
    for(region:listRegion){
        if(region.size<minSize && region.size!=0)
            listSmall.push_back(region);
    }
    sort(listSmall);
    
    // Merge de la liste des petites régions
    while(listSmall.size!=0){
        mergeRegion(currentRegion, voisin0);
        update(listSmall, voisin0);
    }
}

void mergeRegion(int idRegion1, int idRegion2){
    //mise à jour de la nouvelle région
    for(int i = 0; i<listId.size; i++)
        if(listId[i]==idRegion1) listId[i]=idRegion2;
    Region r1 = listRegion[idRegion1];
    Region r2 = listRegion[idRegion2];
    listRegion[idRegion2].avg = (r2.avg * r2.size + r1.moy * r1.size)/(r1.size + r2.size);
    listRegion[idRegion2].size += r1.size;

    //mise à jour de la liste d'adjacence de la nouvelle région
    adjacence[idRegion2].delete(idRegion1);
    adjacence[idRegion1].delete(idRegion2);
    adjacence[idRegion2].fusion(adjacence[idRegion1]);
    sort(adjacence[idRegion2]);
}

void update(vector<region> &list, int idRegionMerge){
    list.erase(0);

    if(idRegionMerge in list){
        sort(list);
    }
}