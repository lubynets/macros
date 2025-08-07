std::pair<float, float> evalArmenterosPodolanskiCascade(aod::CascData const& casc)
{
  const float pxv0 = casc.pxpos() + casc.pxneg();
  const float pyv0 = casc.pypos() + casc.pyneg();
  const float pzv0 = casc.pzpos() + casc.pzneg();

  const float pxbach = casc.pxbach();
  const float pybach = casc.pybach();
  const float pzbach = casc.pzbach();

  const double momTot = RecoDecay::p2(pxv0 + pxbach, pyv0 + pybach, pzv0 + pzbach);

  const float lQlv0 = RecoDecay::dotProd(std::array{pxv0, pyv0, pzv0}, std::array{pxv0 + pxbach, pyv0 + pybach, pzv0 + pzbach}) / momTot;
  const float lQlbach = RecoDecay::dotProd(std::array{pxbach, pybach, pzbach}, std::array{pxv0 + pxbach, pyv0 + pybach, pzv0 + pzbach}) / momTot;
  const float dp = RecoDecay::dotProd(std::array{pxbach, pybach, pzbach}, std::array{pxv0 + pxbach, pyv0 + pybach, pzv0 + pzbach});

  float alpha = (lQlv0 - lQlbach) / (lQlv0 + lQlbach);
  const float qtarm = std::sqrt(RecoDecay::p2(pxneg, pyneg, pzneg) - dp * dp / momTot / momTot);

  if(casc.sign() > 0) { // TODO check
    alpha = -alpha;
  }

  return std::make_pair(alpha, qtarm);
}

Configurable<float> cutAPOmegaUp1{"cutAPOmegaUp1", 0.25, "cutAPOmegaUp1"};
Configurable<float> cutAPOmegaUp2{"cutAPOmegaUp2", 0.358, "cutAPOmegaUp2"};
Configurable<float> cutAPOmegaUp3{"cutAPOmegaUp3", 0.35, "cutAPOmegaUp3"};
Configurable<float> cutAPOmegaDown1{"cutAPOmegaDown1", 0.15, "cutAPOmegaDown1"};
Configurable<float> cutAPOmegaDown2{"cutAPOmegaDown2", 0.358, "cutAPOmegaDown2"};
Configurable<float> cutAPOmegaDown3{"cutAPOmegaDown3", 0.16, "cutAPOmegaDown3"};
Configurable<float> cutAlphaOmegaHigh{"cutAlphaOmegaHigh", 0.358, "cutAlphaOmegaHigh"};
Configurable<float> cutAlphaOmegaLow{"cutAlphaOmegaLow", 0., "cutAlphaOmegaLow"};
Configurable<float> cutMassOmegaHigh{"cutMassOmegaHigh", 1.677, "cutMassOmegaHigh"};
Configurable<float> cutMassOmegaLow{"cutMassOmegaLow", 1.667, "cutMassOmegaLow"};

if(cutAPOmegaUp3 > cutAPOmegaUp2) {
  LOGP(error, "cutAPOmegaUp3 must be less than cutAPOmegaUp2");
}

registry.add("hCascAPplotSelected", "hCascAPplotSelected", HistType::kTH2F, {{200, -1.0f, +1.0f}, {250, 0.0f, 0.25f}});
registry.add("hMassOmega", "hMassOmega", HistType::kTH2F, {{900, 0.0f, 90.0f}, {100, 1.62f, 1.72f}});
registry.add("hMassAntiOmega", "hMassAntiOmega", HistType::kTH2F, {{900, 0.0f, 90.0f}, {100, 1.62f, 1.72f}});

int checkCascade(float alpha, float qt)
{
  const bool isAlphaPos = alpha > 0;
  alpha = std::fabs(alpha);

  const float qUp = cutAPOmegaUp1 * std::sqrt(std::abs(1.0f - ((alpha - cutAPOmegaUp2) * (alpha - cutAPOmegaUp2)) / (cutAPOmegaUp3 * cutAPOmegaUp3)));
  const float qDown = cutAPOmegaDown1 * std::sqrt(std::abs(1.0f - ((alpha - cutAPOmegaDown2) * (alpha - cutAPOmegaDown2)) / (cutAPOmegaDown3 * cutAPOmegaDown3)));

  if(alpha < cutAlphaOmegaLow || alpha > cutAlphaOmegaHigh || qt < qDown || qt > qUp) {
    return kUndef;
  }

  if(isAlphaPos) {
    return kOmega;
  } else {
    return kAntiOmega;
  }
}

for (const auto& Casc : Cascs) {
  if (fillhisto) {
    registry.fill(HIST("hCascCandidate"), 1);
  }

  // topological selections

  if (fillhisto) {
    registry.fill(HIST("hCascCandidate"), 2);
  }

  const float Cascradius = Casc.cascradius();

  const float mOmega = Casc.mOmega();

  auto [alpha, qt] = evalArmenterosPodolanskiCascade(Casc);

  const int cascid = checkCascade(alpha, qt);
  if (cascid < 0) {
    continue;
  }
  if (fillhisto) {
    registry.fill(HIST("hCascAPplotSelected"), alpha, qt);
  }

  auto storeCascAddID = [&](auto gix, auto id) {
    if (produceCascID.value) {
      cascpidmap[gix] = id;
    }
  };

  if(cascid == kOmega) {
    if (fillhisto) {
      registry.fill(HIST("hMassOmega"), Cascradius, mOmega);
    }
    if(cutMassOmegaLow < mOmega && mOmega < cutMassOmegaHigh) {
      pidmap[Casc.bachelorId()] |= (uint8_t(1) << kOmega);
      storeCascAddID(Casc.globalIndex(), kOmega);
    }
  } else if (
    if (fillhisto) {
      registry.fill(HIST("hMassAntiOmega"), Cascradius, mOmega);
    }
    if(cutMassOmegaLow < mOmega && mOmega < cutMassOmegaHigh) {
      pidmap[Casc.bachelorId()] |= (uint8_t(1) << kAntiOmega);
      storeCascAddID(Casc.globalIndex(), kAntiOmega);
    }
  )



} // end of Casc loop




































