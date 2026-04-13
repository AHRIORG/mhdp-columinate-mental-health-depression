const STORAGE_KEY = "columinate-cmi-prototype-v6";
const ALL_DOMAINS = "__all_domains__";
const WORDING_MODES = {
  improved: "AI improved",
  original: "Original"
};

const domainPriority = [
  "Poverty and inequality",
  "Education and employment",
  "Housing and living conditions",
  "Early childhood development",
  "Social support and exclusion",
  "Physical environment",
  "Basic services",
  "Violence and safety",
  "Religion and spirituality"
];

function defaultPlaceForDomain(domain) {
  const lookup = {
    "Poverty and inequality": "Household and livelihood setting",
    "Education and employment": "School, work, and household setting",
    "Housing and living conditions": "Household and neighborhood setting",
    "Early childhood development": "Early learning and care setting",
    "Social support and exclusion": "Family, peer, and community setting",
    "Physical environment": "Home and environmental setting",
    "Basic services": "Clinic, infrastructure, and service setting",
    "Violence and safety": "Home, tavern, and community setting",
    "Religion and spirituality": "Faith and community setting"
  };
  return lookup[domain] || "Community setting";
}

const baseSampleSpecs = [
  { id: "poverty-food", domain: "Poverty and inequality", title: "Food insecurity pathway", pathway: "Poverty -> inability to afford food -> stress and worry -> mental health decline -> Poor mental health." },
  { id: "subsistence-loss-chaos", domain: "Poverty and inequality", title: "Subsistence farming loss pathway", pathway: "Reliance on subsistence farming -> loss of income due to theft or natural hazards -> household chaos -> anger and domestic violence -> negative mental health outcomes -> Poor mental health." },
  { id: "poverty-home-frustration", domain: "Poverty and inequality", title: "Poverty and living conditions pathway", pathway: "Poverty -> inability to fix or improve living conditions independently -> anxiety and frustration -> Poor mental health." },
  { id: "hope-future", domain: "Poverty and inequality", title: "Hope for the future pathway", pathway: "Hope for the future -> encouragement to persist -> reduced hopelessness -> positive outlook -> improved mental health -> Positive mental health." },

  { id: "school-fencing-risk", domain: "Education and employment", title: "Unsafe school environment pathway", pathway: "Lack of school fencing and monitoring -> truancy and time in unsafe spaces -> substance use and sexual activity -> emotional strain and disrupted education -> Poor mental health." },
  { id: "school-punishment-humiliation", domain: "Education and employment", title: "Physical punishment at school pathway", pathway: "Physical punishment at school -> humiliation -> suppressed emotions -> anger and resentment toward authority -> Poor mental health." },
  { id: "mining-income", domain: "Education and employment", title: "Mining employment pathway", pathway: "Employment through mining -> financial security -> improved mental health -> Positive mental health." },
  { id: "mine-closure-stress", domain: "Education and employment", title: "Mine closure pathway", pathway: "Mine closure -> loss of income -> reliance on grants -> communal stress and anxiety -> Poor mental health." },
  { id: "financial-stress-substance", domain: "Education and employment", title: "Financial stress and substance use pathway", pathway: "Financial stress -> emotional strain -> substance use -> family conflict -> emotional distress -> Poor mental health." },
  { id: "unemployment-self-esteem", domain: "Education and employment", title: "Unemployment pathway", pathway: "Protracted unemployment -> eroded self-esteem -> stress -> Poor mental health." },
  { id: "grant-dependence-self-esteem", domain: "Education and employment", title: "Grant dependence pathway", pathway: "Reliance on government grant -> eroded self-esteem -> Poor mental health." },

  { id: "water-hygiene-shame", domain: "Housing and living conditions", title: "Water and hygiene pathway", pathway: "Lack of running water -> inability to maintain hygiene -> bullying and shame -> emotional distress -> Poor mental health." },
  { id: "poor-housing-strain", domain: "Housing and living conditions", title: "Poor housing pathway", pathway: "Poor condition of houses -> discomfort at home -> mental strain -> Poor mental health." },

  { id: "early-learning-self-esteem", domain: "Early childhood development", title: "Early learning performance pathway", pathway: "Access to early learning centers -> better performance in higher grades -> increased self-esteem and reduced stress -> Positive mental health." },
  { id: "early-learning-safety", domain: "Early childhood development", title: "Early learning safety pathway", pathway: "Access to early learning centers -> safety -> less exposure to trauma -> Positive mental health." },

  { id: "migration-loneliness", domain: "Social support and exclusion", title: "Migration and loneliness pathway", pathway: "Migration for work -> family separation -> loneliness and insecurity -> increased stress and anxiety -> Poor mental health." },
  { id: "church-guidance-peace", domain: "Social support and exclusion", title: "Church guidance pathway", pathway: "Church counselling and spiritual guidance -> sense of hope and peace -> Positive mental health." },
  { id: "sports-engagement-protection", domain: "Social support and exclusion", title: "Sports engagement pathway", pathway: "Access to sports grounds -> social engagement through sport -> reduced substance abuse -> improved mental health -> Positive mental health." },
  { id: "bad-peers-conflict", domain: "Social support and exclusion", title: "Peer conflict pathway", pathway: "Association with bad peers -> parental disapproval -> conflict -> emotional strain -> Poor mental health." },
  { id: "exclusion-decision-making", domain: "Social support and exclusion", title: "Family decision-making exclusion pathway", pathway: "Feeling excluded from family decision-making -> feelings of marginalisation and powerlessness -> Poor mental health." },
  { id: "community-social-capital", domain: "Social support and exclusion", title: "Community social capital pathway", pathway: "Strong family and community social capital -> access to shared resources such as land, gardens, and employment opportunities -> material security and social support -> improved well-being -> Positive mental health." },
  { id: "low-connectedness-exclusion", domain: "Social support and exclusion", title: "Low social connectedness pathway", pathway: "Low social connectedness and weak community cohesion -> exclusion from shared resources and fear of stigma -> reduced help-seeking and exposure to harm such as abuse or conflict -> stress and emotional distress -> Poor mental health." },

  { id: "pollution-fear", domain: "Physical environment", title: "Pollution and worry pathway", pathway: "Air and water pollution -> ill health -> fear and worry -> Poor mental health." },
  { id: "mine-land-grief", domain: "Physical environment", title: "Mine land damage pathway", pathway: "Mines -> land damage -> grief and stress -> mental health deterioration -> Poor mental health." },
  { id: "dumping-site-avoidance", domain: "Physical environment", title: "Dumping site pathway", pathway: "Litter from dumping sites -> perceived danger -> avoidance of nature -> reduced relaxation and well-being -> Poor mental health." },
  { id: "subsistence-agriculture-relief", domain: "Physical environment", title: "Subsistence agriculture pathway", pathway: "Subsistence agriculture -> lower financial pressure -> Positive mental health." },
  { id: "extreme-weather-loss", domain: "Physical environment", title: "Extreme weather pathway", pathway: "Extreme weather events -> soil erosion and crop destruction -> loss of food and income -> mental distress -> Poor mental health." },
  { id: "water-scarcity-conflict", domain: "Physical environment", title: "Water scarcity pathway", pathway: "Water scarcity -> dual use of water such as drinking water for cattle versus watering crops -> internal conflict -> emotional distress -> Poor mental health." },

  { id: "infrastructure-danger", domain: "Basic services", title: "Infrastructure danger pathway", pathway: "Neglected public infrastructure -> exposure to physical danger -> persistent worry and compensatory costs -> Poor mental health." },
  { id: "infrastructure-missed-opportunities", domain: "Basic services", title: "Infrastructure access pathway", pathway: "Neglected public infrastructure -> inability to use services -> missed opportunities -> frustration -> Poor mental health." },
  { id: "political-will-distress", domain: "Basic services", title: "Low political will pathway", pathway: "Low political will to improve basic services -> unmet basic needs -> sense of injustice -> emotional distress -> Poor mental health." },
  { id: "health-awareness-relief", domain: "Basic services", title: "Health awareness pathway", pathway: "Improved health awareness -> lower stigma -> easier access to treatment -> emotional relief and improved physical health -> Positive mental health." },
  { id: "clinic-queues", domain: "Basic services", title: "Clinic queue pathway", pathway: "Long queues at clinics -> delays in or avoidance of care -> increased anxiety and discouragement -> Poor mental health." },
  { id: "stigmatizing-treatment", domain: "Basic services", title: "Stigmatising treatment pathway", pathway: "Stigmatising treatment from healthcare workers -> humiliation -> avoidance of treatment -> worsening health and mental strain -> Poor mental health." },
  { id: "no-nearby-hospital", domain: "Basic services", title: "Hospital distance pathway", pathway: "Lack of nearby hospital -> lengthy travel time and transport poverty -> hinders urgent and intensive medical care -> sense of distance from needed care -> Poor mental health." },
  { id: "domestic-violence-withdrawal", domain: "Basic services", title: "Domestic violence case withdrawal pathway", pathway: "Domestic violence -> women withdraw cases due to financial dependence -> ongoing distress -> Poor mental health." },
  { id: "flooded-bridge-disruption", domain: "Basic services", title: "Flooded bridge disruption pathway", pathway: "Neglected public infrastructure -> unable to cross bridges during floods -> unable to attend school or work -> negatively impacts education and income -> Poor mental health." },

  { id: "taverns-risk-violence", domain: "Violence and safety", title: "Taverns and unsafe spaces pathway", pathway: "Abandoned buildings and taverns -> drug and alcohol use -> increased violence, crime, and risky behaviours -> emotional distress and insecurity -> Poor mental health." },
  { id: "substance-domestic-violence", domain: "Violence and safety", title: "Substance use and domestic violence pathway", pathway: "Drug and alcohol use -> domestic violence -> Poor mental health." },
  { id: "substance-debt-poverty", domain: "Violence and safety", title: "Substance use and debt pathway", pathway: "Drug and alcohol use -> misuse of income -> debt and poverty -> Poor mental health." },

  { id: "church-attendance-support", domain: "Religion and spirituality", title: "Church attendance pathway", pathway: "Church attendance and pastoral counselling -> support, mediation, and guidance -> Positive mental health." },
  { id: "prayer-rituals-relief", domain: "Religion and spirituality", title: "Prayer and rituals pathway", pathway: "Prayer and rituals -> spiritual connection and perceived ancestral support -> relief and hope -> Positive mental health." },
  { id: "witchcraft-fear", domain: "Religion and spirituality", title: "Fear of witchcraft pathway", pathway: "Fear of witchcraft -> worry and distrust -> mental distress -> Poor mental health." },
  { id: "dual-belief-conflict", domain: "Religion and spirituality", title: "Dual belief systems pathway", pathway: "Dual belief systems -> conflicting norms -> internal conflict -> emotional strain -> Poor mental health." },
  { id: "church-rejection-shame", domain: "Religion and spirituality", title: "Church rejection pathway", pathway: "Rejection by church when engaged in stigmatised behaviour -> loss of status -> feeling shunned -> emotional distress -> Poor mental health." }
];

const originalPathwaysById = {
  "poverty-food": "Poverty -> inability to afford food -> stress and worry -> mental health decline -> Poor mental health",
  "subsistence-loss-chaos": "Reliance on subsistence farming -> loss of income due to theft or natural hazards -> household chaos -> anger and domestic violence -> negative mental health outcomes -> Poor mental health",
  "poverty-home-frustration": "Poverty -> inability to fix or improve living conditions independently -> anxiety and frustration -> Poor mental health",
  "hope-future": "Hope for the future -> encouragement to persist -> reduced hopelessness -> positive outlook -> improved mental health -> Positive mental health",
  "school-fencing-risk": "Lack of school fencing and monitoring -> truancy and time in unsafe spaces -> substance use and sexual activity -> emotional strain and disrupted education -> Poor mental health",
  "school-punishment-humiliation": "Physical punishment at school -> humiliation -> suppressed emotions -> anger and resentment toward authority -> Poor mental health",
  "mining-income": "Employment through mining -> financial security -> improved mental health -> Positive mental health",
  "mine-closure-stress": "Mine closure -> loss of income -> reliance on grants -> communal stress and anxiety -> Poor mental health",
  "financial-stress-substance": "Financial stress -> emotional strain -> substance use -> family conflict -> emotional distress -> Poor mental health",
  "unemployment-self-esteem": "Protracted unemployment -> eroded self-esteem -> stress -> Poor mental health",
  "grant-dependence-self-esteem": "Reliance on government grant -> eroded self-esteem -> Poor mental health",
  "water-hygiene-shame": "Lack of running water -> inability to maintain hygiene -> bullying and shame -> emotional distress -> Poor mental health",
  "poor-housing-strain": "Poor condition of houses -> discomfort at home -> mental strain -> Poor mental health",
  "early-learning-self-esteem": "Access to early learning centers -> better performance in higher grades -> increased self-esteem and reduced stress -> Positive mental health",
  "early-learning-safety": "Access to early learning centers -> safety -> less exposure to trauma -> Positive mental health",
  "migration-loneliness": "Migration for work -> family separation -> loneliness and insecurity -> increased stress and anxiety -> Poor mental health",
  "church-guidance-peace": "Church counselling and spiritual guidance -> sense of hope and peace -> Positive mental health",
  "sports-engagement-protection": "Access to sports grounds -> social engagement through sport -> reduced substance abuse -> improved mental health -> Positive mental health",
  "bad-peers-conflict": "Association with bad peers -> parental disapproval -> conflict -> emotional strain -> Poor mental health",
  "exclusion-decision-making": "Feeling excluded from family decision-making -> feelings of marginalisation and powerlessness -> Poor mental health",
  "community-social-capital": "Strong family and community social capital -> access to shared resources such as land, gardens, and employment opportunities -> material security and social support -> improved well-being -> Positive mental health",
  "low-connectedness-exclusion": "Low social connectedness and weak community cohesion -> exclusion from shared resources and fear of stigma -> reduced help-seeking and exposure to harm such as abuse or conflict -> stress and emotional distress -> Poor mental health",
  "pollution-fear": "Air and water pollution -> ill health -> fear and worry -> Poor mental health",
  "mine-land-grief": "Mines -> land damage -> grief and stress -> mental health deterioration -> Poor mental health",
  "dumping-site-avoidance": "Litter from dumping sites -> perceived danger -> avoidance of nature -> reduced relaxation and well-being -> Poor mental health",
  "subsistence-agriculture-relief": "Subsistence agriculture -> lower financial pressure -> Positive mental health",
  "extreme-weather-loss": "Extreme weather events -> soil erosion and crop destruction -> loss of food and income -> mental distress -> Poor mental health",
  "water-scarcity-conflict": "Water scarcity -> dual use of water such as drinking water for cattle versus watering crops -> internal conflict -> emotional distress -> Poor mental health",
  "infrastructure-danger": "Neglected public infrastructure -> exposure to physical danger -> persistent worry and compensatory costs -> Poor mental health",
  "infrastructure-missed-opportunities": "Neglected public infrastructure -> inability to use services -> missed opportunities -> frustration -> Poor mental health",
  "political-will-distress": "Low political will to improve basic services -> unmet basic needs -> sense of injustice -> emotional distress -> Poor mental health",
  "health-awareness-relief": "Improved health awareness -> lower stigma -> easier access to treatment -> emotional relief and improved physical health -> Positive mental health",
  "clinic-queues": "Long queues at clinics -> delays in or avoidance of care -> increased anxiety and discouragement -> Poor mental health",
  "stigmatizing-treatment": "Stigmatising treatment from healthcare workers -> humiliation -> avoidance of treatment -> worsening health and mental strain -> Poor mental health",
  "no-nearby-hospital": "Lack of nearby hospital -> lengthy travel time and transport poverty -> hinders urgent and intensive medical care -> sense of distance from needed care -> Poor mental health",
  "domestic-violence-withdrawal": "Domestic violence -> women withdraw cases due to financial dependence -> ongoing distress -> Poor mental health",
  "flooded-bridge-disruption": "Neglected public infrastructure -> unable to cross bridges during floods -> unable to attend school or work -> negatively impacts education and income -> Poor mental health",
  "taverns-risk-violence": "Abandoned buildings and taverns -> drug and alcohol use -> increased violence, crime, and risky behaviours -> emotional distress and insecurity -> Poor mental health",
  "substance-domestic-violence": "Drug and alcohol use -> domestic violence -> Poor mental health",
  "substance-debt-poverty": "Drug and alcohol use -> misuse of income -> debt and poverty -> Poor mental health",
  "church-attendance-support": "Church attendance and pastoral counselling -> support, mediation, and guidance -> Positive mental health",
  "prayer-rituals-relief": "Prayer and rituals -> spiritual connection and perceived ancestral support -> relief and hope -> Positive mental health",
  "witchcraft-fear": "Fear of witchcraft -> worry and distrust -> mental distress -> Poor mental health",
  "dual-belief-conflict": "Dual belief systems -> conflicting norms -> internal conflict -> emotional strain -> Poor mental health",
  "church-rejection-shame": "Rejection by church when engaged in stigmatised behaviour -> loss of status -> feeling shunned -> emotional distress -> Poor mental health"
};

function defaultNodes(polarity) {
  return [
    { type: "driver", label: "New driver" },
    { type: "bridge", label: "Bridge" },
    { type: "outcome", label: outcomeLabel(polarity) }
  ];
}

function createWordingVersion(pathway, story, polarity) {
  const normalizedPathway = String(pathway || "").replace(/\s+/g, " ").trim();
  const labels = normalizedPathway.replace(/\.$/, "").split(/\s*->\s*/).filter(Boolean);
  const nodes = labels.length > 0 ? classifyNodeLabels(labels, polarity) : defaultNodes(polarity);
  return {
    story: String(story || normalizedPathway || "").trim(),
    nodes
  };
}

function harmonizeNodeLabel(label, polarity) {
  const replacementRules = [
    [/\bcounselling\b/gi, "counseling"],
    [/\bstigmatising\b/gi, "stigmatizing"],
    [/\bmarginalisation\b/gi, "marginalization"],
    [/\bbehaviours\b/gi, "behaviors"],
    [/\bwell-being\b/gi, "wellbeing"],
    [/\bmental health decline\b/gi, "psychological distress"],
    [/\bmental health deterioration\b/gi, "psychological distress"],
    [/\bnegative mental health outcomes\b/gi, "psychological distress"],
    [/\bill health\b/gi, "physical ill health"]
  ];

  let cleaned = String(label || "").replace(/\s+/g, " ").trim();
  replacementRules.forEach(([pattern, value]) => {
    cleaned = cleaned.replace(pattern, value);
  });

  if (/^poor mental health$/i.test(cleaned)) {
    return "Higher risk of poor mental health";
  }
  if (/^positive mental health$/i.test(cleaned)) {
    return "Stronger positive mental health";
  }

  return cleaned;
}

function harmonizePathwayText(pathway, polarity) {
  const labels = String(pathway || "").replace(/\s+/g, " ").trim().replace(/\.$/, "").split(/\s*->\s*/).filter(Boolean);
  if (labels.length === 0) return "";

  const harmonized = labels.map((label) => harmonizeNodeLabel(label, polarity));
  harmonized[harmonized.length - 1] = polarity === "positive"
    ? "Stronger positive mental health"
    : "Higher risk of poor mental health";
  return harmonized.join(" -> ");
}

function pathwaySample(spec, index) {
  const improvedSourcePathway = String(spec.pathway || "").replace(/\s+/g, " ").trim();
  const improvedPathway = harmonizePathwayText(improvedSourcePathway, /positive mental health\.?$/i.test(improvedSourcePathway) ? "positive" : "negative");
  const polarity = /positive mental health\.?$/i.test(improvedPathway) ? "positive" : "negative";
  const improvedLabels = improvedPathway.replace(/\.$/, "").split(/\s*->\s*/).filter(Boolean);
  const originalPathway = String(originalPathwaysById[spec.id] || improvedPathway).replace(/\s+/g, " ").trim();

  return {
    id: spec.id || slugify(`inventory-${index + 1}-${improvedLabels[0]}`),
    domain: spec.domain,
    title: spec.title || `${improvedLabels[0]} pathway`,
    place: spec.place || defaultPlaceForDomain(spec.domain),
    polarity,
    reviewStatus: "approved",
    wordingVersions: {
      improved: createWordingVersion(improvedPathway, spec.story || improvedPathway, polarity),
      original: createWordingVersion(originalPathway, spec.originalStory || originalPathway, polarity)
    }
  };
}

const baseSamples = baseSampleSpecs.map(pathwaySample);
const samples = JSON.parse(JSON.stringify(baseSamples));

const state = {
  activeDomain: ALL_DOMAINS,
  activeId: samples[0].id,
  activeTab: "overview",
  query: "",
  polarityFilter: "all",
  wordingMode: "improved",
  draft: null,
  isNew: false
};

const bannerTitle = document.getElementById("bannerTitle");
const currentStatus = document.getElementById("currentStatus");
const chainSentence = document.getElementById("chainSentence");
const domainList = document.getElementById("domainList");
const savedList = document.getElementById("savedList");
const savedMeta = document.getElementById("savedMeta");
const chainTitle = document.getElementById("chainTitle");
const chainSummary = document.getElementById("chainSummary");
const chainCanvas = document.getElementById("chainCanvas");
const titleInput = document.getElementById("titleInput");
const domainSelect = document.getElementById("domainSelect");
const polaritySelect = document.getElementById("polaritySelect");
const placeInput = document.getElementById("placeInput");
const storyInput = document.getElementById("storyInput");
const nodeList = document.getElementById("nodeList");
const reviewQueue = document.getElementById("reviewQueue");
const queueSummary = document.getElementById("queueSummary");
const searchInput = document.getElementById("searchInput");
const overviewTotal = document.getElementById("overviewTotal");
const overviewDraft = document.getElementById("overviewDraft");
const overviewApproved = document.getElementById("overviewApproved");
const overviewCurrent = document.getElementById("overviewCurrent");
const overviewSentence = document.getElementById("overviewSentence");
const dataStatus = document.getElementById("dataStatus");
const importFile = document.getElementById("importFile");
const wordingModeHelp = document.getElementById("wordingModeHelp");
const builderModeCopy = document.getElementById("builderModeCopy");

function slugify(value) {
  return String(value || "untitled")
    .toLowerCase()
    .replace(/[^a-z0-9]+/g, "-")
    .replace(/^-+|-+$/g, "") || `chain-${Date.now()}`;
}

function deepClone(value) {
  return JSON.parse(JSON.stringify(value));
}

function normalizeWordingMode(value) {
  return value === "original" ? "original" : "improved";
}

function wordingModeLabel(value) {
  return WORDING_MODES[normalizeWordingMode(value)];
}

function outcomeLabel(polarity) {
  return polarity === "positive" ? "Positive mental health" : "Poor mental health";
}

function normalizePolarity(value) {
  return String(value || "").toLowerCase().includes("positive") ? "positive" : "negative";
}

function normalizeReviewStatus(value) {
  const candidate = String(value || "draft").toLowerCase();
  if (["draft", "reviewed", "approved"].includes(candidate)) return candidate;
  return "draft";
}

function statusLabel(status) {
  return status.charAt(0).toUpperCase() + status.slice(1);
}

function sentenceFromNodes(nodes) {
  return nodes.map((node) => String(node.label || "").trim()).filter(Boolean).join(" -> ");
}

function normalizeNodes(nodes, polarity) {
  const clean = Array.isArray(nodes)
    ? nodes.map((node) => ({
      type: node?.type,
      label: String(node?.label || node?.name || "").trim()
    })).filter((node) => node.label)
    : [];

  if (clean.length === 0) {
    return defaultNodes(polarity);
  }

  return clean.map((node, index) => {
    const inferredType = ["driver", "bridge", "outcome"].includes(node.type)
      ? node.type
      : (index === 0 ? "driver" : index === clean.length - 1 ? "outcome" : "bridge");
    return {
      type: inferredType,
      label: inferredType === "outcome" ? outcomeLabel(polarity) : node.label
    };
  });
}

function normalizeWordingVersion(rawVersion, fallbackVersion, polarity) {
  const fallback = fallbackVersion || createWordingVersion("", "", polarity);
  const sourceNodes = Array.isArray(rawVersion?.nodes) && rawVersion.nodes.length > 0
    ? rawVersion.nodes
    : fallback.nodes;
  const nodes = normalizeNodes(sourceNodes, polarity);
  const fallbackStory = fallback.story || sentenceFromNodes(nodes);
  return {
    story: String(rawVersion?.story ?? fallbackStory ?? "").trim(),
    nodes
  };
}

function activeVersionFor(chain, wordingMode = state.wordingMode) {
  const mode = normalizeWordingMode(wordingMode);
  if (!chain.wordingVersions) {
    chain.wordingVersions = {
      improved: createWordingVersion(chain.story || sentenceFromNodes(chain.nodes || defaultNodes(chain.polarity)), chain.story, chain.polarity),
      original: createWordingVersion(chain.story || sentenceFromNodes(chain.nodes || defaultNodes(chain.polarity)), chain.story, chain.polarity)
    };
  }
  if (!chain.wordingVersions[mode]) {
    const fallbackMode = mode === "original" ? "improved" : "original";
    chain.wordingVersions[mode] = deepClone(chain.wordingVersions[fallbackMode] || createWordingVersion("", "", chain.polarity));
  }
  return chain.wordingVersions[mode];
}

function prepareWordingVersions(chain) {
  const fallbackSentence = sentenceFromNodes(chain.nodes || defaultNodes(chain.polarity));
  const fallbackVersion = createWordingVersion(fallbackSentence, chain.story || fallbackSentence, chain.polarity);
  chain.wordingVersions = {
    improved: normalizeWordingVersion(chain.wordingVersions?.improved, fallbackVersion, chain.polarity),
    original: normalizeWordingVersion(chain.wordingVersions?.original, fallbackVersion, chain.polarity)
  };
}

function syncOutcomeLabels(chain) {
  prepareWordingVersions(chain);
  Object.values(chain.wordingVersions).forEach((version) => {
    version.nodes = normalizeNodes(version.nodes, chain.polarity);
  });
}

function sentenceFor(chain, wordingMode = state.wordingMode) {
  const labels = activeVersionFor(chain, wordingMode).nodes.map((node) => node.label.trim()).filter(Boolean);
  return labels.length > 0 ? labels.join(" -> ") : "Start building the chain by adding a driver.";
}

function summaryFor(chain) {
  return `${chain.domain} pathway ending in ${outcomeLabel(chain.polarity).toLowerCase()}.`;
}

function classifyNodeLabels(labels, polarity) {
  const clean = labels.map((label) => String(label || "").trim()).filter(Boolean);
  return clean.map((label, index) => {
    if (index === 0) return { type: "driver", label };
    if (index === clean.length - 1) return { type: "outcome", label: outcomeLabel(polarity) };
    return { type: "bridge", label };
  });
}

function normalizeSample(raw, index = 0) {
  const polarity = normalizePolarity(raw.polarity || raw.outcome || raw.outcome_polarity);
  const reviewStatus = normalizeReviewStatus(raw.reviewStatus || raw.review_status);
  let nodes = [];

  if (Array.isArray(raw.nodes) && raw.nodes.length > 0) {
    nodes = normalizeNodes(raw.nodes, polarity);
  } else if (raw.pathway_sentence) {
    nodes = classifyNodeLabels(String(raw.pathway_sentence).split(/\s*->\s*/), polarity);
  }

  if (nodes.length === 0) {
    nodes = defaultNodes(polarity);
  }

  const legacySentence = String(
    raw.pathway_sentence ||
    raw.ai_improved_pathway_sentence ||
    raw.improved_pathway_sentence ||
    raw.original_pathway_sentence ||
    sentenceFromNodes(nodes)
  ).trim();
  const legacyVersion = createWordingVersion(
    legacySentence,
    raw.story || raw.story_note || raw.evidence || raw.ai_improved_story_note || raw.original_story_note || legacySentence,
    polarity
  );

  const originalRawVersion = raw.wordingVersions?.original || raw.versions?.original || (
    raw.original_pathway_sentence || raw.original_story_note
      ? createWordingVersion(raw.original_pathway_sentence || legacySentence, raw.original_story_note || raw.story || raw.story_note, polarity)
      : null
  );
  const improvedRawVersion = raw.wordingVersions?.improved || raw.versions?.improved || (
    raw.ai_improved_pathway_sentence || raw.improved_pathway_sentence || raw.ai_improved_story_note
      ? createWordingVersion(
        raw.ai_improved_pathway_sentence || raw.improved_pathway_sentence || legacySentence,
        raw.ai_improved_story_note || raw.story || raw.story_note,
        polarity
      )
      : null
  );

  return {
    id: slugify(raw.id || raw.chain_id || raw.title || raw.chain_title || `imported-${index + 1}`),
    domain: raw.domain || domainPriority[0],
    title: raw.title || raw.chain_title || `Imported pathway ${index + 1}`,
    place: raw.place || raw.location || "",
    polarity,
    reviewStatus,
    wordingVersions: {
      improved: normalizeWordingVersion(improvedRawVersion, legacyVersion, polarity),
      original: normalizeWordingVersion(originalRawVersion, legacyVersion, polarity)
    }
  };
}

function matchesQuery(sample, query) {
  if (!query) return true;
  prepareWordingVersions(sample);
  const wordingText = Object.values(sample.wordingVersions).flatMap((version) => [
    version.story,
    ...version.nodes.map((node) => node.label)
  ]);
  const haystack = [
    sample.title,
    sample.domain,
    sample.place,
    ...wordingText
  ].join(" ").toLowerCase();
  return haystack.includes(query.toLowerCase());
}

function filteredSamples() {
  return samples.filter((sample) => {
    const polarityOk = state.polarityFilter === "all" || sample.polarity === state.polarityFilter;
    return polarityOk && matchesQuery(sample, state.query);
  });
}

function samplesForDomain(domain) {
  if (domain === ALL_DOMAINS) {
    return filteredSamples();
  }
  return filteredSamples().filter((sample) => sample.domain === domain);
}

function availableDomains() {
  return domainPriority.filter((domain) => samplesForDomain(domain).length > 0);
}

function setStatus(message) {
  dataStatus.textContent = message;
}

function saveBrowserState() {
  try {
    window.localStorage.setItem(STORAGE_KEY, JSON.stringify({
      samples,
      state: {
        activeDomain: state.activeDomain,
        activeId: state.activeId,
        activeTab: state.activeTab,
        query: state.query,
        polarityFilter: state.polarityFilter,
        wordingMode: state.wordingMode,
        draft: state.draft,
        isNew: state.isNew
      }
    }));
  } catch (_error) {
    // no-op
  }
}

function loadBrowserState() {
  try {
    const raw = window.localStorage.getItem(STORAGE_KEY);
    if (!raw) return;
    const parsed = JSON.parse(raw);
    if (Array.isArray(parsed?.samples)) {
      samples.splice(0, samples.length, ...parsed.samples.map(normalizeSample));
    }
    if (parsed?.state && typeof parsed.state === "object") {
      Object.assign(state, parsed.state);
    }
    state.wordingMode = normalizeWordingMode(state.wordingMode);
  } catch (_error) {
    // no-op
  }
}

function resetBrowserState() {
  try {
    window.localStorage.removeItem(STORAGE_KEY);
  } catch (_error) {
    // no-op
  }
  window.location.reload();
}

function ensureActiveDomain() {
  const domains = availableDomains();
  if (filteredSamples().length === 0) {
    state.activeDomain = ALL_DOMAINS;
    return;
  }
  if (state.activeDomain === ALL_DOMAINS) {
    return;
  }
  if (!domains.includes(state.activeDomain)) {
    state.activeDomain = ALL_DOMAINS;
  }
}

function ensureActiveChain() {
  const current = samples.find((sample) => sample.id === state.activeId);
  if (current && matchesQuery(current, state.query) &&
      (state.polarityFilter === "all" || current.polarity === state.polarityFilter)) {
    return;
  }
  const domainItems = samplesForDomain(state.activeDomain);
  if (domainItems.length > 0) {
    state.activeId = domainItems[0].id;
    loadDraft(state.activeId);
  }
}

function populateDomainSelect() {
  domainSelect.innerHTML = "";
  domainPriority.forEach((domain) => {
    const option = document.createElement("option");
    option.value = domain;
    option.textContent = domain;
    domainSelect.appendChild(option);
  });
}

function loadDraft(id) {
  const found = samples.find((sample) => sample.id === id) || samples[0];
  const preserveAllDomains = state.activeDomain === ALL_DOMAINS;
  state.activeId = found.id;
  state.activeDomain = preserveAllDomains ? ALL_DOMAINS : found.domain;
  state.draft = normalizeSample(deepClone(found));
  state.isNew = false;
}

function createBlankChain() {
  const domain = state.activeDomain && state.activeDomain !== ALL_DOMAINS
    ? state.activeDomain
    : domainPriority[0];
  state.draft = {
    id: `draft-${Date.now()}`,
    domain,
    title: "New pathway",
    place: "",
    polarity: "negative",
    reviewStatus: "draft",
    wordingVersions: {
      improved: createWordingVersion("", "", "negative"),
      original: createWordingVersion("", "", "negative")
    }
  };
  state.activeId = state.draft.id;
  state.isNew = true;
  state.activeTab = "builder";
}

function syncDetailInputs() {
  const activeVersion = activeVersionFor(state.draft);
  titleInput.value = state.draft.title;
  domainSelect.value = state.draft.domain;
  polaritySelect.value = state.draft.polarity;
  placeInput.value = state.draft.place;
  storyInput.value = activeVersion.story;
  searchInput.value = state.query;
}

function renderTabs() {
  const validTabs = new Set(["overview", "library", "builder", "review"]);
  if (!validTabs.has(state.activeTab)) {
    state.activeTab = "overview";
  }

  document.querySelectorAll(".tab-btn").forEach((button) => {
    const active = button.dataset.tab === state.activeTab;
    button.classList.toggle("active", active);
    button.setAttribute("aria-selected", active ? "true" : "false");
  });

  document.querySelectorAll(".tab-panel").forEach((panel) => {
    const active = panel.id === `panel-${state.activeTab}`;
    panel.hidden = !active;
    panel.classList.toggle("active", active);
  });
}

function renderOverview() {
  overviewTotal.textContent = String(samples.length);
  overviewDraft.textContent = String(samples.filter((sample) => sample.reviewStatus === "draft").length);
  overviewApproved.textContent = String(samples.filter((sample) => sample.reviewStatus === "approved").length);
  overviewCurrent.textContent = state.draft.title;
  overviewSentence.textContent = sentenceFor(state.draft);
}

function renderDomains() {
  ensureActiveDomain();
  domainList.innerHTML = "";
  const domains = availableDomains();
  const totalVisible = filteredSamples().length;

  if (totalVisible === 0) {
    domainList.innerHTML = '<div class="empty-note">No pathways match the current search.</div>';
    return;
  }

  const allButton = document.createElement("button");
  allButton.type = "button";
  allButton.className = `domain-btn${state.activeDomain === ALL_DOMAINS ? " active" : ""}`;
  allButton.innerHTML = `<strong>All domains</strong><small>${totalVisible} matching pathways</small>`;
  allButton.addEventListener("click", () => {
    state.activeDomain = ALL_DOMAINS;
    renderAll();
  });
  domainList.appendChild(allButton);

  domains.forEach((domain) => {
    const button = document.createElement("button");
    button.type = "button";
    button.className = `domain-btn${domain === state.activeDomain ? " active" : ""}`;
    button.innerHTML = `<strong>${domain}</strong><small>${samplesForDomain(domain).length} matching pathways</small>`;
    button.addEventListener("click", () => {
      state.activeDomain = domain;
      const first = samplesForDomain(domain)[0];
      if (first) loadDraft(first.id);
      renderAll();
    });
    domainList.appendChild(button);
  });
}

function renderSaved() {
  const domainItems = samplesForDomain(state.activeDomain);
  savedList.innerHTML = "";
  const scopeText = state.activeDomain === ALL_DOMAINS
    ? `${domainItems.length} pathway${domainItems.length === 1 ? "" : "s"} across all domains in this view.`
    : `${domainItems.length} pathway${domainItems.length === 1 ? "" : "s"} in this view.`;
  savedMeta.textContent = `${scopeText} Showing ${wordingModeLabel(state.wordingMode).toLowerCase()} sentences.`;

  if (domainItems.length === 0) {
    savedList.innerHTML = '<div class="empty-note">No saved pathways match the current filters.</div>';
    return;
  }

  domainItems.forEach((sample) => {
    const button = document.createElement("button");
    button.type = "button";
    button.className = `saved-btn${sample.id === state.activeId ? " active" : ""}`;
    button.innerHTML = `
      <strong>${sample.title}</strong>
      <small>${summaryFor(sample)}</small>
      <span class="saved-sentence">${sentenceFor(sample)}</span>
      <small>${sample.reviewStatus}</small>
    `;
    button.addEventListener("click", () => {
      loadDraft(sample.id);
      state.activeTab = "builder";
      renderAll();
    });
    savedList.appendChild(button);
  });
}

function renderCanvas() {
  chainCanvas.innerHTML = "";
  const activeVersion = activeVersionFor(state.draft);
  activeVersion.nodes.forEach((node, index) => {
    const card = document.createElement("div");
    const outcomeClass = node.type === "outcome" ? ` outcome ${state.draft.polarity}` : "";
    card.className = `node-card ${node.type}${outcomeClass}`;
    card.innerHTML = `<strong>${node.type}</strong>${node.label || "Untitled node"}`;
    chainCanvas.appendChild(card);

    if (index < activeVersion.nodes.length - 1) {
      const arrow = document.createElement("div");
      arrow.className = "arrow";
      arrow.textContent = "→";
      chainCanvas.appendChild(arrow);
    }
  });
}

function renderNodeEditor() {
  nodeList.innerHTML = "";
  const activeVersion = activeVersionFor(state.draft);
  activeVersion.nodes.forEach((node, index) => {
    const row = document.createElement("div");
    row.className = "node-row";

    const typeSelect = document.createElement("select");
    ["driver", "bridge", "outcome"].forEach((type) => {
      const option = document.createElement("option");
      option.value = type;
      option.textContent = type.charAt(0).toUpperCase() + type.slice(1);
      option.selected = node.type === type;
      typeSelect.appendChild(option);
    });
    typeSelect.addEventListener("change", (event) => {
      activeVersion.nodes[index].type = event.target.value;
      if (event.target.value === "outcome") {
        activeVersion.nodes[index].label = outcomeLabel(state.draft.polarity);
      }
      renderWorkspace();
    });

    const labelInput = document.createElement("input");
    labelInput.type = "text";
    labelInput.value = node.label;
    labelInput.placeholder = "Node label";
    labelInput.addEventListener("input", (event) => {
      activeVersion.nodes[index].label = event.target.value;
      renderWorkspaceMeta();
    });

    const actions = document.createElement("div");
    actions.className = "node-actions";

    const moveUp = document.createElement("button");
    moveUp.type = "button";
    moveUp.className = "ghost-btn";
    moveUp.textContent = "↑";
    moveUp.disabled = index === 0;
    moveUp.addEventListener("click", () => {
      const temp = activeVersion.nodes[index - 1];
      activeVersion.nodes[index - 1] = activeVersion.nodes[index];
      activeVersion.nodes[index] = temp;
      renderWorkspace();
    });

    const moveDown = document.createElement("button");
    moveDown.type = "button";
    moveDown.className = "ghost-btn";
    moveDown.textContent = "↓";
    moveDown.disabled = index === activeVersion.nodes.length - 1;
    moveDown.addEventListener("click", () => {
      const temp = activeVersion.nodes[index + 1];
      activeVersion.nodes[index + 1] = activeVersion.nodes[index];
      activeVersion.nodes[index] = temp;
      renderWorkspace();
    });

    const remove = document.createElement("button");
    remove.type = "button";
    remove.className = "ghost-btn";
    remove.textContent = "Remove";
    remove.disabled = activeVersion.nodes.length <= 2;
    remove.addEventListener("click", () => {
      activeVersion.nodes.splice(index, 1);
      renderWorkspace();
    });

    actions.append(moveUp, moveDown, remove);
    row.append(typeSelect, labelInput, actions);
    nodeList.appendChild(row);
  });
}

function renderWorkspaceMeta(syncInputs = false) {
  const wordingLabel = wordingModeLabel(state.wordingMode);
  bannerTitle.textContent = state.draft.title;
  currentStatus.textContent = statusLabel(state.draft.reviewStatus);
  currentStatus.className = `state-badge ${state.draft.reviewStatus}`;
  chainSentence.textContent = sentenceFor(state.draft);
  chainTitle.textContent = state.draft.title;
  chainSummary.textContent = summaryFor(state.draft);
  wordingModeHelp.textContent = `Switch between the cleaned AI-improved wording and the original pathway inventory wording. Currently viewing: ${wordingLabel}.`;
  builderModeCopy.textContent = `Editing the ${wordingLabel.toLowerCase()} wording for this pathway.`;
  if (syncInputs) syncDetailInputs();
  renderCanvas();
  renderOverview();
  saveBrowserState();
}

function renderWorkspace() {
  renderWorkspaceMeta(true);
  renderNodeEditor();
}

function renderQueue() {
  const counts = ["draft", "reviewed", "approved"].map((status) => ({
    status,
    count: samples.filter((sample) => sample.reviewStatus === status).length
  }));

  queueSummary.innerHTML = "";
  counts.forEach((entry) => {
    const card = document.createElement("div");
    card.className = "queue-card";
    card.innerHTML = `<span>${statusLabel(entry.status)}</span><strong>${entry.count}</strong>`;
    queueSummary.appendChild(card);
  });

  reviewQueue.innerHTML = "";
  samples
    .slice()
    .sort((a, b) => ({ draft: 0, reviewed: 1, approved: 2 }[a.reviewStatus] - { draft: 0, reviewed: 1, approved: 2 }[b.reviewStatus]))
    .forEach((sample) => {
      const button = document.createElement("button");
      button.type = "button";
      button.className = `queue-item${sample.id === state.activeId ? " active" : ""}`;
      button.innerHTML = `
        <span class="queue-title">${sample.title}</span>
        <span class="queue-meta">${sample.domain}</span>
        <span class="queue-state ${sample.reviewStatus}">${statusLabel(sample.reviewStatus)}</span>
      `;
      button.addEventListener("click", () => {
        loadDraft(sample.id);
        state.activeTab = "builder";
        renderAll();
      });
      reviewQueue.appendChild(button);
    });
}

function renderAll() {
  ensureActiveDomain();
  ensureActiveChain();
  document.getElementById("polarityFilters").querySelectorAll("button").forEach((chip) => {
    chip.classList.toggle("active", chip.dataset.filter === state.polarityFilter);
  });
  document.querySelectorAll("[data-wording]").forEach((button) => {
    button.classList.toggle("active", button.dataset.wording === state.wordingMode);
  });
  renderTabs();
  renderDomains();
  renderSaved();
  renderWorkspace();
  renderQueue();
  saveBrowserState();
}

function persistDraft(nextStatus = null) {
  const activeVersion = activeVersionFor(state.draft);
  state.draft.title = titleInput.value.trim() || "Untitled pathway";
  state.draft.domain = domainSelect.value;
  state.draft.polarity = polaritySelect.value;
  state.draft.place = placeInput.value.trim();
  activeVersion.story = storyInput.value.trim();
  syncOutcomeLabels(state.draft);
  prepareWordingVersions(state.draft);
  Object.values(state.draft.wordingVersions).forEach((version) => {
    version.story = String(version.story || sentenceFromNodes(version.nodes)).trim();
    version.nodes = version.nodes.map((node) => {
      if (node.type === "outcome") return { type: "outcome", label: outcomeLabel(state.draft.polarity) };
      return { type: node.type, label: node.label.trim() || `${node.type} step` };
    });
  });
  if (nextStatus) state.draft.reviewStatus = nextStatus;

  const sampleIndex = samples.findIndex((sample) => sample.id === state.draft.id);
  const saved = deepClone(state.draft);
  if (sampleIndex >= 0) {
    samples[sampleIndex] = saved;
  } else {
    saved.id = slugify(saved.title);
    state.draft.id = saved.id;
    samples.unshift(saved);
  }

  state.activeId = state.draft.id;
  state.activeDomain = state.draft.domain;
  state.isNew = false;
  setStatus(`Saved ${saved.title}.`);
}

function setDraftFromInputs() {
  const activeVersion = activeVersionFor(state.draft);
  state.draft.title = titleInput.value;
  state.draft.domain = domainSelect.value;
  state.draft.polarity = polaritySelect.value;
  state.draft.place = placeInput.value;
  activeVersion.story = storyInput.value;
  syncOutcomeLabels(state.draft);
}

function csvEscape(value) {
  return `"${String(value ?? "").replace(/"/g, '""')}"`;
}

function buildEdgeRows() {
  const rows = [["chain_id", "chain_title", "domain", "polarity", "review_status", "wording_mode", "step_order", "source_label", "target_label"]];
  samples.forEach((sample) => {
    const activeVersion = activeVersionFor(sample);
    for (let i = 0; i < activeVersion.nodes.length - 1; i += 1) {
      rows.push([
        sample.id,
        sample.title,
        sample.domain,
        sample.polarity,
        sample.reviewStatus,
        wordingModeLabel(state.wordingMode),
        i + 1,
        activeVersion.nodes[i].label,
        activeVersion.nodes[i + 1].label
      ]);
    }
  });
  return rows.map((row) => row.map(csvEscape).join(",")).join("\n");
}

function buildSummaryRows() {
  const rows = [[
    "chain_id",
    "title",
    "domain",
    "polarity",
    "review_status",
    "place",
    "selected_wording",
    "display_pathway_sentence",
    "original_pathway_sentence",
    "ai_improved_pathway_sentence",
    "display_story_note",
    "original_story_note",
    "ai_improved_story_note"
  ]];
  samples.forEach((sample) => {
    const originalVersion = activeVersionFor(sample, "original");
    const improvedVersion = activeVersionFor(sample, "improved");
    const displayVersion = activeVersionFor(sample);
    rows.push([
      sample.id,
      sample.title,
      sample.domain,
      sample.polarity,
      sample.reviewStatus,
      sample.place,
      wordingModeLabel(state.wordingMode),
      sentenceFromNodes(displayVersion.nodes),
      sentenceFromNodes(originalVersion.nodes),
      sentenceFromNodes(improvedVersion.nodes),
      displayVersion.story,
      originalVersion.story,
      improvedVersion.story
    ]);
  });
  return rows.map((row) => row.map(csvEscape).join(",")).join("\n");
}

function buildStateJson() {
  return JSON.stringify({ samples, state }, null, 2);
}

function downloadFile(filename, content, mimeType) {
  const blob = new Blob([content], { type: mimeType });
  const url = URL.createObjectURL(blob);
  const link = document.createElement("a");
  link.href = url;
  link.download = filename;
  document.body.appendChild(link);
  link.click();
  link.remove();
  URL.revokeObjectURL(url);
}

function parseCsv(text) {
  const rows = [];
  let row = [];
  let value = "";
  let inQuotes = false;

  for (let i = 0; i < text.length; i += 1) {
    const char = text[i];
    const next = text[i + 1];

    if (char === '"') {
      if (inQuotes && next === '"') {
        value += '"';
        i += 1;
      } else {
        inQuotes = !inQuotes;
      }
    } else if (char === "," && !inQuotes) {
      row.push(value);
      value = "";
    } else if ((char === "\n" || char === "\r") && !inQuotes) {
      if (char === "\r" && next === "\n") i += 1;
      row.push(value);
      if (row.some((cell) => cell !== "")) rows.push(row);
      row = [];
      value = "";
    } else {
      value += char;
    }
  }

  row.push(value);
  if (row.some((cell) => cell !== "")) rows.push(row);
  return rows;
}

function buildSamplesFromSummaryCsv(rows) {
  const headers = rows[0].map((header) => header.trim().toLowerCase());
  return rows.slice(1).map((row, index) => {
    const record = {};
    headers.forEach((header, headerIndex) => {
      record[header] = row[headerIndex] || "";
    });
    return normalizeSample(record, index);
  });
}

function buildSamplesFromEdgeCsv(rows) {
  const headers = rows[0].map((header) => header.trim().toLowerCase());
  const map = new Map();
  rows.slice(1).forEach((row) => {
    const record = {};
    headers.forEach((header, headerIndex) => {
      record[header] = row[headerIndex] || "";
    });
    const key = record.chain_id || record.chain_title || `edge-${map.size + 1}`;
    if (!map.has(key)) map.set(key, []);
    map.get(key).push(record);
  });

  return Array.from(map.values()).map((records, index) => {
    records.sort((a, b) => Number(a.step_order || 0) - Number(b.step_order || 0));
    const labels = [];
    records.forEach((record, recordIndex) => {
      if (recordIndex === 0) labels.push(record.source_label);
      labels.push(record.target_label);
    });
    const first = records[0] || {};
    return normalizeSample({
      chain_id: first.chain_id,
      chain_title: first.chain_title,
      domain: first.domain,
      polarity: first.polarity,
      review_status: first.review_status,
      pathway_sentence: labels.join(" -> ")
    }, index);
  });
}

function applyImportedSamples(importedSamples, importedState = null) {
  if (!importedSamples.length) throw new Error("No valid pathways were found in the selected file.");

  const replace = window.confirm("Click OK to replace the current browser data with the imported file. Click Cancel to append imported pathways instead.");
  if (replace) {
    samples.splice(0, samples.length, ...importedSamples);
  } else {
    importedSamples.forEach((sample) => {
      const existingIndex = samples.findIndex((entry) => entry.id === sample.id);
      if (existingIndex >= 0) samples[existingIndex] = sample;
      else samples.unshift(sample);
    });
  }

  if (importedState && replace) {
    Object.assign(state, importedState);
  }

  state.wordingMode = normalizeWordingMode(state.wordingMode);
  state.activeId = importedSamples[0].id;
  state.activeDomain = importedSamples[0].domain;
  state.activeTab = "builder";
  loadDraft(state.activeId);
  setStatus(`Imported ${importedSamples.length} pathway${importedSamples.length === 1 ? "" : "s"}.`);
  renderAll();
}

function handleImportText(name, text) {
  const lowerName = name.toLowerCase();
  if (lowerName.endsWith(".json")) {
    const parsed = JSON.parse(text);
    if (Array.isArray(parsed)) {
      applyImportedSamples(parsed.map(normalizeSample));
      return;
    }
    if (Array.isArray(parsed.samples)) {
      applyImportedSamples(parsed.samples.map(normalizeSample), parsed.state || null);
      return;
    }
    if (Array.isArray(parsed.chains)) {
      applyImportedSamples(parsed.chains.map(normalizeSample));
      return;
    }
    throw new Error("Unsupported JSON structure.");
  }

  const rows = parseCsv(text);
  if (rows.length < 2) throw new Error("CSV file does not contain any data rows.");
  const headers = rows[0].map((header) => header.trim().toLowerCase());
  if (
    headers.includes("pathway_sentence") ||
    headers.includes("display_pathway_sentence") ||
    headers.includes("ai_improved_pathway_sentence") ||
    headers.includes("original_pathway_sentence")
  ) {
    applyImportedSamples(buildSamplesFromSummaryCsv(rows));
    return;
  }
  if (headers.includes("source_label") && headers.includes("target_label")) {
    applyImportedSamples(buildSamplesFromEdgeCsv(rows));
    return;
  }
  throw new Error("Unsupported CSV structure. Use summary CSV or edge list CSV exported by this app.");
}

function bindFieldListeners() {
  document.querySelectorAll(".tab-btn").forEach((button) => {
    button.addEventListener("click", () => {
      state.activeTab = button.dataset.tab;
      renderTabs();
      saveBrowserState();
    });
  });

  document.querySelectorAll("[data-go-tab]").forEach((button) => {
    button.addEventListener("click", () => {
      state.activeTab = button.dataset.goTab;
      renderTabs();
      saveBrowserState();
    });
  });

  document.querySelectorAll("[data-wording]").forEach((button) => {
    button.addEventListener("click", () => {
      state.wordingMode = normalizeWordingMode(button.dataset.wording);
      renderAll();
    });
  });

  [titleInput, domainSelect, polaritySelect, placeInput, storyInput].forEach((element) => {
    element.addEventListener("input", () => {
      setDraftFromInputs();
      if (element === polaritySelect) renderWorkspace();
      else renderWorkspaceMeta();
    });
    element.addEventListener("change", () => {
      setDraftFromInputs();
      if (element === polaritySelect) renderWorkspace();
      else renderWorkspaceMeta();
    });
  });

  searchInput.addEventListener("input", (event) => {
    state.query = event.target.value.trim();
    renderAll();
  });

  document.getElementById("polarityFilters").querySelectorAll("button").forEach((button) => {
    button.addEventListener("click", () => {
      state.polarityFilter = button.dataset.filter;
      document.getElementById("polarityFilters").querySelectorAll("button").forEach((chip) => {
        chip.classList.toggle("active", chip === button);
      });
      renderAll();
    });
  });

  document.getElementById("btnNewChain").addEventListener("click", () => {
    createBlankChain();
    renderAll();
  });

  document.getElementById("btnAddDriver").addEventListener("click", () => {
    activeVersionFor(state.draft).nodes.splice(0, 0, { type: "driver", label: "New driver" });
    renderWorkspace();
  });

  document.getElementById("btnAddBridge").addEventListener("click", () => {
    const activeVersion = activeVersionFor(state.draft);
    const outcomeIndex = activeVersion.nodes.findIndex((node) => node.type === "outcome");
    const insertAt = outcomeIndex >= 0 ? outcomeIndex : activeVersion.nodes.length;
    activeVersion.nodes.splice(insertAt, 0, { type: "bridge", label: "New bridge" });
    renderWorkspace();
  });

  document.getElementById("btnAddOutcome").addEventListener("click", () => {
    const activeVersion = activeVersionFor(state.draft);
    const existingOutcome = activeVersion.nodes.findIndex((node) => node.type === "outcome");
    if (existingOutcome >= 0) activeVersion.nodes[existingOutcome].label = outcomeLabel(state.draft.polarity);
    else activeVersion.nodes.push({ type: "outcome", label: outcomeLabel(state.draft.polarity) });
    renderWorkspace();
  });

  document.getElementById("btnSaveDraft").addEventListener("click", () => {
    persistDraft("draft");
    loadDraft(state.activeId);
    renderAll();
  });

  document.getElementById("btnMarkReviewed").addEventListener("click", () => {
    persistDraft("reviewed");
    loadDraft(state.activeId);
    renderAll();
  });

  document.getElementById("btnApprove").addEventListener("click", () => {
    persistDraft("approved");
    loadDraft(state.activeId);
    renderAll();
  });

  document.getElementById("btnReset").addEventListener("click", () => {
    if (samples.some((sample) => sample.id === state.activeId)) loadDraft(state.activeId);
    else createBlankChain();
    renderAll();
  });

  document.getElementById("btnExportState").addEventListener("click", () => {
    downloadFile("cmi_state.json", buildStateJson(), "application/json;charset=utf-8");
    setStatus("Downloaded JSON state export.");
  });

  document.getElementById("btnExportEdges").addEventListener("click", () => {
    downloadFile("cmi_edge_list.csv", buildEdgeRows(), "text/csv;charset=utf-8");
    setStatus("Downloaded edge list CSV.");
  });

  document.getElementById("btnExportSummaries").addEventListener("click", () => {
    downloadFile("cmi_chain_summaries.csv", buildSummaryRows(), "text/csv;charset=utf-8");
    setStatus("Downloaded chain summary CSV.");
  });

  document.getElementById("btnImportData").addEventListener("click", () => {
    importFile.value = "";
    importFile.click();
  });

  importFile.addEventListener("change", async (event) => {
    const file = event.target.files?.[0];
    if (!file) return;
    try {
      const text = await file.text();
      handleImportText(file.name, text);
    } catch (error) {
      setStatus(error.message || "Import failed.");
    }
  });

  document.getElementById("btnResetBrowser").addEventListener("click", () => {
    resetBrowserState();
  });
}

loadBrowserState();
populateDomainSelect();
if (state.draft && state.draft.id && samples.some((sample) => sample.id === state.draft.id)) {
  state.draft = normalizeSample(state.draft);
} else {
  loadDraft(state.activeId);
}
bindFieldListeners();
setStatus("Supported imports: JSON state exports, summary CSV, or edge list CSV. Summary exports now include both original and AI-improved wording.");
renderAll();
